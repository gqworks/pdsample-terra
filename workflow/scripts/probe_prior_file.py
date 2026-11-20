#!/usr/bin/env python3
"""Decide whether a Terra prior path is a real non-empty file.

WDL must not coerce a placeholder/missing gs:// String to File and then size()
it (Cromwell FileNotFoundException at disk_gb eval). This script takes a String
path, never fails the process on a missing object, and writes true/false.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request
import zlib

BLANK = frozenset({"", ".", "null", "none", "na", "n/a", "undefined", "nil", "none."})
INTROME_MARKER = "introme_score"
VCF_PREFIXES = ("##fileformat=vcf", "#chrom", "##fileformat=vcf4")


def normalize_path(path: str) -> str:
    text = (path or "").strip().strip('"').strip("'")
    if text.startswith("file://"):
        text = text[7:]
    return text


def is_blank(path: str) -> bool:
    return normalize_path(path).lower() in BLANK


def gcs_parts(path: str) -> tuple[str, str] | tuple[None, None]:
    if not path.startswith("gs://"):
        return None, None
    rest = path[5:]
    bucket, slash, obj = rest.partition("/")
    if not bucket or not slash or not obj:
        return None, None
    return bucket, obj


def metadata_token() -> str | None:
    req = urllib.request.Request(
        "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
        headers={"Metadata-Flavor": "Google"},
    )
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            payload = json.loads(resp.read().decode("utf-8"))
        token = payload.get("access_token")
        return str(token) if token else None
    except Exception:
        return None


def gcs_request(path: str, token: str | None, method: str = "HEAD", extra_headers: dict | None = None):
    bucket, obj = gcs_parts(path)
    if not bucket:
        raise ValueError("not a gs:// object")
    quoted = urllib.parse.quote(obj, safe="/")
    url = f"https://storage.googleapis.com/{bucket}/{quoted}"
    headers = dict(extra_headers or {})
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = urllib.request.Request(url, method=method, headers=headers)
    return urllib.request.urlopen(req, timeout=30)


def gcs_size(path: str, token: str | None) -> int | None:
    try:
        with gcs_request(path, token, method="HEAD") as resp:
            length = resp.headers.get("Content-Length")
            return int(length) if length is not None else None
    except urllib.error.HTTPError as err:
        if err.code == 404:
            return None
        return None
    except Exception:
        return None


def gcs_peek(path: str, token: str | None, nbytes: int = 65536) -> bytes:
    try:
        with gcs_request(path, token, method="GET", extra_headers={"Range": f"bytes=0-{nbytes - 1}"}) as resp:
            return resp.read(nbytes)
    except Exception:
        return b""


def local_size(path: str) -> int | None:
    try:
        if not os.path.isfile(path):
            return None
        return os.path.getsize(path)
    except OSError:
        return None


def local_peek(path: str, nbytes: int = 65536) -> bytes:
    try:
        with open(path, "rb") as handle:
            return handle.read(nbytes)
    except OSError:
        return b""


def decode_head(data: bytes) -> str:
    if not data:
        return ""
    raw = data
    if data.startswith(b"\x1f\x8b"):
        try:
            raw = gzip.decompress(data)
        except Exception:
            try:
                decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
                raw = decoder.decompress(data)
            except Exception:
                return ""
    return raw.decode("utf-8", "replace")


def looks_introme(data: bytes) -> bool:
    text = decode_head(data)
    for line in text.splitlines():
        if line.startswith("#") and "introme_score" not in line:
            continue
        if not line.strip():
            continue
        return INTROME_MARKER in line.split("\t") or INTROME_MARKER in line
    return False


def looks_vcf(data: bytes) -> bool:
    text = decode_head(data)
    for line in text.splitlines():
        if not line.strip():
            continue
        low = line.lstrip("\ufeff").lower()
        return low.startswith(VCF_PREFIXES) or low.startswith("##")
    return False


def usable(path: str, kind: str = "any", min_bytes: int = 1) -> bool:
    """True only when the object exists, is non-empty, and matches kind."""
    kind = (kind or "any").strip().lower()
    path = normalize_path(path)
    if is_blank(path):
        return False
    if path.startswith("gs://"):
        token = metadata_token()
        size = gcs_size(path, token)
        head = gcs_peek(path, token) if size else b""
    else:
        size = local_size(path)
        head = local_peek(path) if size else b""
    if size is None or size < min_bytes:
        return False
    if kind in {"introme", "introme_tsv"}:
        return looks_introme(head)
    if kind in {"vcf", "vep", "vep_vcf"}:
        return looks_vcf(head)
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path", default="", help="Terra String prior path (not a WDL File)")
    parser.add_argument("--kind", default="any", choices=("any", "introme", "introme_tsv", "vcf", "vep", "vep_vcf"))
    parser.add_argument("--min-bytes", type=int, default=1)
    parser.add_argument("--out", default="usable.txt")
    args = parser.parse_args(argv)
    result = usable(args.path, kind=args.kind, min_bytes=args.min_bytes)
    with open(args.out, "w", encoding="utf-8") as handle:
        handle.write("true\n" if result else "false\n")
    sys.stdout.write("true\n" if result else "false\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
