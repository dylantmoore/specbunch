"""Slow estimator checks are opt-in locally and enabled in full CI."""
import os
import pytest

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False,
                     help="run expensive estimator checks")

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: expensive estimator check")

def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow") or os.environ.get("SPECBUNCH_RUN_SLOW") == "1":
        return
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(pytest.mark.skip(reason="use --runslow for the full suite"))
