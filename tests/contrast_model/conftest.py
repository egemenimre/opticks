# opticks Models and analysis tools for optical system engineering
#
# Copyright (C) Egemen Imre
#
# Licensed under GNU GPL v3.0. See LICENSE.md for more info.
"""
Shared pytest fixtures for the contrast-model test suite.

Provides ``pushbroom_yaml`` — a session-scoped factory fixture that
resolves a yaml filename via the pushbroom-data search paths, eliminating
the repeated path-constant boilerplate from individual test classes.
"""

from pathlib import Path

import pytest

from tests import process_paths

_PUSHBR_DIR = Path("sat_pushbroom_data")
_PUSHBR_ALT_DIR = Path("tests", "imaging_model", "sat_pushbroom_data")


@pytest.fixture(scope="session")
def pushbroom_yaml():
    """Return a callable that resolves a yaml filename via the pushbroom
    search paths (current dir, then the imaging_model fallback)."""

    def _resolve(filename: str) -> Path:
        return process_paths(Path(filename), _PUSHBR_DIR, _PUSHBR_ALT_DIR)

    return _resolve
