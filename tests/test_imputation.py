import numpy as np
import pandas as pd
import pytest

from LIMBR.imputation import imputable


def test_init_loads_data(imputation_file):
    obj = imputable(imputation_file, missingness=0.5)
    assert obj.data is not None
    assert obj.miss == 0.5
    assert obj.pats == {}
    assert obj.notdone is True


def test_init_custom_neighbors(imputation_file):
    obj = imputable(imputation_file, missingness=0.5, neighbors=7)
    assert obj.NN == 7


def test_deduplicate_creates_multiindex(imputation_file):
    obj = imputable(imputation_file, missingness=0.5)
    obj.deduplicate()
    assert isinstance(obj.data.index, pd.MultiIndex)
    assert obj.data.index.names == ["Peptide", "Protein"]


def test_deduplicate_no_duplicate_peptides(imputation_file):
    obj = imputable(imputation_file, missingness=0.5)
    obj.deduplicate()
    counts = obj.data.groupby(level="Peptide").size()
    assert (counts == 1).all()


def test_drop_missing_removes_high_missing_rows(imputation_file):
    obj = imputable(imputation_file, missingness=0.0)  # keep only complete rows
    obj.deduplicate()
    obj.drop_missing()
    # only the 50 complete rows should remain
    assert obj.data.isnull().sum().sum() == 0
    assert len(obj.data) == 50


def test_drop_missing_high_threshold_keeps_all(imputation_file):
    obj = imputable(imputation_file, missingness=1.0)
    obj.deduplicate()
    n_before = len(obj.data)
    obj.drop_missing()
    assert len(obj.data) == n_before


def test_impute_fills_all_nan(imputation_file, tmp_path):
    out = str(tmp_path / "imputed.txt")
    obj = imputable(imputation_file, missingness=0.5, neighbors=5)
    obj.deduplicate()
    obj.drop_missing()
    obj.impute(out)
    result = pd.read_csv(out, sep="\t", index_col=[0, 1])
    assert result.isnull().sum().sum() == 0


def test_impute_preserves_complete_row_count(imputation_file, tmp_path):
    out = str(tmp_path / "imputed.txt")
    obj = imputable(imputation_file, missingness=0.5, neighbors=5)
    obj.deduplicate()
    obj.drop_missing()
    n_rows = len(obj.data)
    obj.impute(out)
    result = pd.read_csv(out, sep="\t", index_col=[0, 1])
    assert len(result) == n_rows


def test_impute_data_pipeline(imputation_file, tmp_path):
    out = str(tmp_path / "pipeline.txt")
    obj = imputable(imputation_file, missingness=0.5, neighbors=5)
    obj.impute_data(out)
    result = pd.read_csv(out, sep="\t", index_col=[0, 1])
    assert result.isnull().sum().sum() == 0
    assert len(result) > 0
