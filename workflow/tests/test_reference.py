"""CRAM header versus reference dictionary."""

from check_reference import compare_headers


def test_reference_md5_mismatch_and_chr_prefix():
    ok = compare_headers(
        "@SQ\tSN:chr1\tLN:10\tM5:abc\n",
        "@SQ\tSN:chr1\tLN:10\tM5:abc\n@SQ\tSN:chr2\tLN:10\tM5:def\n",
    )
    assert ok["ok"] is True

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
