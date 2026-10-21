"""End-to-end tests of the hoshi Docker image.

Runs real commands inside the container, the way a user would, and checks the
files they leave on the host:

- ``hoshi report-medical`` Savont/EMU (+ metadata) -> HTML and PDF
- ``hoshi convert`` Savont -> species count table

Opt-in: ``--docker`` builds a fresh image first, ``--docker-image TAG`` reuses
one (see tests/conftest.py).
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from hoshi import cli

REPO_ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = REPO_ROOT / "test_data"
EMU_REL = "emu_output/test_ind/sample01/sample01_rel-abundance.tsv"
SAVONT_REL = "savont_output/test_ind/savont-out-sample01"
METADATA_REL = "medical/metadata_only.json"

pytestmark = pytest.mark.docker


def _is_pdf(path: Path) -> bool:
    return path.is_file() and path.read_bytes().startswith(b"%PDF")


def _docker_hoshi(image: str, out_dir: Path, *hoshi_args: str) -> subprocess.CompletedProcess[str]:
    """Run `hoshi <args>` in the image; test_data is /in (ro), out_dir is /out."""
    return subprocess.run(
        [
            "docker", "run", "--rm",
            "--user", f"{os.getuid()}:{os.getgid()}",
            "-v", f"{TEST_DATA}:/in:ro",
            "-v", f"{out_dir}:/out",
            image,
            "hoshi", *hoshi_args,
        ],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )


# ─── report-medical ──────────────────────────────────────────────────


@pytest.fixture(scope="module")
def default_report(docker_image, tmp_path_factory):
    """Savont sample + metadata JSON + PDF, rendered once for the module."""
    out_dir = tmp_path_factory.mktemp("docker_report")
    result = _docker_hoshi(
        docker_image,
        out_dir,
        "report-medical",
        f"/in/{SAVONT_REL}",
        "--input-format", "savont",
        "--metadata", f"/in/{METADATA_REL}",
        "--top", "5",
        "--output", "/out/report_medical.html",
        "--pdf",
    )
    return result, out_dir


def test_report_succeeds(default_report):
    result, _ = default_report
    assert result.returncode == 0, result.stderr


def test_report_renders_metadata_and_organisms(default_report):
    _, out_dir = default_report
    html = (out_dir / "report_medical.html").read_text(encoding="utf-8")
    # From test_data/medical/metadata_only.json.
    assert "16S-2026-000184" in html
    assert "Somsri Chaiyaphum" in html
    assert "Validated bacterial 16S reference database" in html
    # From the Savont classifier output.
    assert "Clostridioides difficile" in html
    assert "PATHOGEN DETECTED" in html


def test_report_writes_pdf(default_report):
    # PDF rendering needs WeasyPrint's native libs, installed by apt in the image.
    _, out_dir = default_report
    assert _is_pdf(out_dir / "report_medical.pdf")


def test_report_outputs_owned_by_host_user(default_report):
    _, out_dir = default_report
    for name in ("report_medical.html", "report_medical.pdf"):
        assert (out_dir / name).stat().st_uid == os.getuid()


def test_report_emu_input_top_limits_organisms(docker_image, tmp_path):
    result = _docker_hoshi(
        docker_image,
        tmp_path,
        "report-medical",
        f"/in/{EMU_REL}",
        "--input-format", "emu",
        "--top", "2",
        "--output", "/out/emu.html",
    )
    assert result.returncode == 0, result.stderr
    html = (tmp_path / "emu.html").read_text(encoding="utf-8")
    assert "Clostridioides difficile" in html
    # --top 2 drops the 3rd-ranked species.
    assert "Lactiplantibacillus plantarum" not in html
    assert not (tmp_path / "emu.pdf").exists()


def test_report_missing_input_exits_1(docker_image, tmp_path):
    # The CLI exit code must come back through `docker run`.
    result = _docker_hoshi(
        docker_image,
        tmp_path,
        "report-medical",
        "/in/nope/missing.tsv",
        "--output", "/out/report.html",
    )
    assert result.returncode == 1
    assert "not found" in result.stderr.lower()
    assert not any(tmp_path.iterdir())


# ─── convert: Savont -> species count table ──────────────────────────


@pytest.fixture(scope="module")
def savont_table(docker_image, tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("docker_convert")
    result = _docker_hoshi(
        docker_image,
        out_dir,
        "convert",
        "--input-format", "savont",
        "--output-format", "table",
        f"/in/{SAVONT_REL}",
        "-o", "/out/species_counts.tsv",
    )
    return result, out_dir / "species_counts.tsv"


def test_convert_savont_table_succeeds(savont_table):
    result, table = savont_table
    assert result.returncode == 0, result.stderr
    assert table.is_file()
    assert table.stat().st_uid == os.getuid()


def test_convert_savont_table_content(savont_table):
    _, table = savont_table
    df = pd.read_csv(table, sep="\t")
    assert {"estimated_count", "sequence_identity", "tax_id", "species"} <= set(df.columns)
    # The Savont sample has 300 reads; tax_id + counts are restored.
    assert df["estimated_count"].sum() == 300
    assert df["tax_id"].notna().all()
    assert "Clostridioides difficile" in set(df["species"])


def test_convert_savont_table_matches_source_tree(savont_table, tmp_path):
    # The image must produce the same table as the code in this checkout;
    # a mismatch means the image is stale or packaged differently.
    _, table = savont_table
    local = tmp_path / "local.tsv"
    code = cli.main(
        [
            "convert",
            "--input-format", "savont",
            "--output-format", "table",
            str(TEST_DATA / SAVONT_REL),
            "-o", str(local),
        ]
    )
    assert code == 0
    assert table.read_text() == local.read_text()
