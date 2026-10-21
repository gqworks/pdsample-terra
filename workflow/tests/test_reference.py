"""BAM header versus reference dictionary."""

from pathlib import Path

from check_reference import compare_headers, main


def test_reference_md5_mismatch_and_chr_prefix():
    ok = compare_headers(
        "@SQ\tSN:chr1\tLN:10\tM5:abc\n",
        "@SQ\tSN:chr1\tLN:10\tM5:abc\n@SQ\tSN:chr2\tLN:10\tM5:def\n",
    )
    assert ok["ok"] is True
    assert ok["bam_contigs"] == 1

    bad = compare_headers(
        "@SQ\tSN:chr1\tLN:10\tM5:abc\n",
        "@SQ\tSN:chr1\tLN:10\tM5:zzz\n",
    )
    assert bad["ok"] is False
    assert bad["md5_mismatch"][0]["contig"] == "chr1"

    prefix = compare_headers(
        "@SQ\tSN:1\tLN:10\tM5:abc\n",
        "@SQ\tSN:chr1\tLN:10\tM5:abc\n",
    )
    assert prefix["ok"] is False
    assert "chr" in prefix["prefix_hint"]


def test_cli_accepts_bam_header_and_cram_header_alias(tmp_path: Path):
    header = tmp_path / "header.txt"
    ref = tmp_path / "ref.dict"
    header.write_text("@SQ\tSN:chr1\tLN:10\tM5:abc\n", encoding="utf-8")
    ref.write_text("@SQ\tSN:chr1\tLN:10\tM5:abc\n", encoding="utf-8")

    bam_out = tmp_path / "bam.json"
    assert (
        main(
            [
                "--bam-header",
                str(header),
                "--reference-dict",
                str(ref),
                "--out",
                str(bam_out),
            ]
        )
        == 0
    )
    assert '"bam_contigs": 1' in bam_out.read_text(encoding="utf-8")

    cram_out = tmp_path / "cram.json"
    assert (
        main(
            [
                "--cram-header",
                str(header),
                "--reference-dict",
                str(ref),
                "--out",
                str(cram_out),
            ]
        )
        == 0
    )
    assert '"bam_contigs": 1' in cram_out.read_text(encoding="utf-8")
