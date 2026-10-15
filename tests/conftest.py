from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA_DIR


@pytest.fixture(scope="session")
def pdb_path_4men() -> Path:
    return DATA_DIR / "raw" / "4MEN.pdb"


@pytest.fixture(scope="session")
def pdb_path_4meo() -> Path:
    return DATA_DIR / "raw" / "4MEO.pdb"
