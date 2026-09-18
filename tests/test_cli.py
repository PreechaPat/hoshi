from hoshi import cli

def test_no_args(capsys):
    exit_code = cli.main([])
    captured = capsys.readouterr()

    assert exit_code == 1
    assert "usage:" in captured.out
