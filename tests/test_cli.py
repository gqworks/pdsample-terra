"""The Assess task calls ``pdvar assess sample.json``. That name must stay a subcommand."""

from typer.testing import CliRunner

from pdvar.cli import app

runner = CliRunner()


def test_assess_remains_a_subcommand():
    listed = runner.invoke(app, ["--help"])
    assert listed.exit_code == 0
    assert "assess" in listed.stdout

    help_assess = runner.invoke(app, ["assess", "--help"])
    assert help_assess.exit_code == 0
    assert "sample_json" in help_assess.stdout.lower() or "SAMPLE_JSON" in help_assess.stdout
    assert "unexpected extra argument" not in help_assess.output.lower()
