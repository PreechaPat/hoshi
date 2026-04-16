from hoshi import cli


def test_version_flag(capsys):
    exit_code = cli.main(["--version"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert captured.err == ""
    assert captured.out.strip() == f"hoshi {cli._resolve_version()}"


def test_short_version_flag(capsys):
    exit_code = cli.main(["-v"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert captured.err == ""
    assert captured.out.strip() == f"hoshi {cli._resolve_version()}"
