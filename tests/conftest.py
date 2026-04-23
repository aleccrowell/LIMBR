import numpy as np
import pandas as pd
import pytest

NTPOINTS = 12
NREPS = 2
NCOLS = NTPOINTS * NREPS  # 24
NROWS = 60


def _col_names():
    return [f"ZT{tp * 2:02d}_{r}" for tp in range(1, NTPOINTS + 1) for r in range(1, NREPS + 1)]


@pytest.fixture
def rnaseq_file(tmp_path):
    np.random.seed(42)
    cols = _col_names()
    df = pd.DataFrame(np.random.randn(NROWS, len(cols)), columns=cols)
    df.index = [f"gene_{i}" for i in range(NROWS)]
    df.index.name = "#"
    path = tmp_path / "rnaseq.txt"
    df.to_csv(path, sep="\t")
    return str(path)


@pytest.fixture
def proteomics_file(tmp_path):
    np.random.seed(42)
    cols = _col_names()
    df = pd.DataFrame(np.random.randn(NROWS, len(cols)), columns=cols)
    df.insert(0, "Protein", [f"prot_{i}" for i in range(NROWS)])
    df.insert(0, "Peptide", [f"pep_{i}" for i in range(NROWS)])
    path = tmp_path / "proteomics.txt"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def imputation_file(tmp_path):
    np.random.seed(42)
    cols = _col_names()
    data = np.random.randn(NROWS, len(cols))
    # last 10 rows get 3 missing values each; first 50 are complete
    for i in range(NROWS - 10, NROWS):
        miss_idx = np.random.choice(len(cols), 3, replace=False)
        data[i, miss_idx] = np.nan
    df = pd.DataFrame(data, columns=cols)
    df.insert(0, "Protein", [f"prot_{i}" for i in range(NROWS)])
    df.insert(0, "Peptide", [f"pep_{i}" for i in range(NROWS)])
    path = tmp_path / "imputation.txt"
    df.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def block_design_file(tmp_path):
    blocks = [0] * (NCOLS // 2) + [1] * (NCOLS // 2)
    path = tmp_path / "blocks.parquet"
    pd.DataFrame({"block": blocks}).to_parquet(str(path))
    return str(path)
