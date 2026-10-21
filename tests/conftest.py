"""Shared pytest setup.

Docker E2E tests (marked ``docker``) are opt-in:

    --docker               build a fresh image from ./Dockerfile, then test it
    --docker-image TAG     test an existing image, no build

Without either flag they are skipped.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
# Own tag, so a test build never replaces hoshi:latest.
E2E_BUILD_TAG = "hoshi:e2e-test"


def pytest_addoption(parser: pytest.Parser) -> None:
    group = parser.getgroup("docker")
    group.addoption(
        "--docker",
        action="store_true",
        help=f"Build {E2E_BUILD_TAG} from ./Dockerfile and run the docker E2E tests.",
    )
    group.addoption(
        "--docker-image",
        metavar="TAG",
        help="Run the docker E2E tests against an existing image (no build).",
    )


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if config.getoption("--docker") or config.getoption("--docker-image"):
        return
    skip = pytest.mark.skip(reason="docker E2E: pass --docker or --docker-image TAG")
    for item in items:
        if "docker" in item.keywords:
            item.add_marker(skip)


@pytest.fixture(scope="session")
def docker_image(pytestconfig: pytest.Config) -> str:
    """Tag of the image under test, building it first when --docker is given."""
    existing = pytestconfig.getoption("--docker-image")
    if existing:
        probe = subprocess.run(
            ["docker", "image", "inspect", existing], capture_output=True, check=False
        )
        if probe.returncode != 0:
            pytest.fail(f"docker image {existing!r} not found; build it or use --docker")
        return existing

    # docker build directly, not scripts/build_docker.sh: that script wipes
    # dist/ and retags hoshi:latest, which a test run should not do.
    build = subprocess.run(
        ["docker", "build", "-t", E2E_BUILD_TAG, "."],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
    )
    if build.returncode != 0:
        log_tail = "\n".join(build.stderr.splitlines()[-40:])
        pytest.fail(f"docker build of {E2E_BUILD_TAG} failed:\n{log_tail}")
    return E2E_BUILD_TAG
