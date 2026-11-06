"""The Assess task calls ``pdvar sample.json`` with no subcommand."""

from typer.testing import CliRunner

from pdvar.cli import app

runner = CliRunner()


def test_sample_json_is_a_direct_argument():
    listed = runner.invoke(app, ["--help"])
    assert listed.exit_code == 0
    assert "sample_json" in listed.stdout.lower() or "SAMPLE_JSON" in listed.stdout
    assert "commands" not in listed.stdout.lower()
