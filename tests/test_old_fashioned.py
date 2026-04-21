import os
import numpy as np
import pandas as pd
import pytest

from LIMBR.old_fashioned import old_fashioned


def test_init_rnaseq(rnaseq_file):
    obj = old_fashioned(rnaseq_file, data_type="r")
    assert obj.data_type == "r"
    assert obj.raw_data is not None
    assert obj.notdone is True


def test_init_proteomics(proteomics_file):
    obj = old_fashioned(proteomics_file, data_type="p")
    assert obj.data_type == "p"
    assert isinstance(obj.raw_data.index, pd.MultiIndex)


def test_pool_normalize_sets_data(rnaseq_file):
    obj = old_fashioned(rnaseq_file, data_type="r")
    obj.pool_normalize()
    assert hasattr(obj, "data")
    assert obj.data.shape == obj.raw_data.shape


def test_pool_normalize_preserves_index(rnaseq_file):
    obj = old_fashioned(rnaseq_file, data_type="r")
    obj.pool_normalize()
    # qnorm sorts rows internally; verify same set of genes present
    assert set(obj.data.index) == set(obj.raw_data.index)


def test_pool_normalize_preserves_columns(rnaseq_file):
    obj = old_fashioned(rnaseq_file, data_type="r")
    obj.pool_normalize()
    assert list(obj.data.columns) == list(obj.raw_data.columns)


def test_normalize_creates_file(rnaseq_file, tmp_path):
    out = str(tmp_path / "normalized.txt")
    obj = old_fashioned(rnaseq_file, data_type="r")
    obj.pool_normalize()
    obj.normalize(out)
    assert os.path.exists(out)


def test_normalize_output_columns(rnaseq_file, tmp_path):
    out = str(tmp_path / "normalized.txt")
    obj = old_fashioned(rnaseq_file, data_type="r")
    obj.pool_normalize()
    obj.normalize(out)
    result = pd.read_csv(out, sep="\t", index_col=0)
    assert list(result.columns) == list(obj.raw_data.columns)


def test_normalize_output_index_name(rnaseq_file, tmp_path):
    out = str(tmp_path / "normalized.txt")
    obj = old_fashioned(rnaseq_file, data_type="r")
    obj.pool_normalize()
    obj.normalize(out)
    result = pd.read_csv(out, sep="\t", index_col=0)
    assert result.index.name == "#"


def test_normalize_output_row_count(rnaseq_file, tmp_path):
    out = str(tmp_path / "normalized.txt")
    obj = old_fashioned(rnaseq_file, data_type="r")
    obj.pool_normalize()
    obj.normalize(out)
    result = pd.read_csv(out, sep="\t", index_col=0)
    assert len(result) == len(obj.raw_data)
