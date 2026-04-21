import os
import numpy as np
import pandas as pd
import pytest

from LIMBR.simulations import simulate


def test_default_shapes():
    sim = simulate(tpoints=12, nrows=50, nreps=2)
    assert len(sim.circ) == 50
    assert sim.sim_miss.shape == (50, 12 * 2)


def test_cols_format():
    sim = simulate(tpoints=12, nrows=50, nreps=2, tpoint_space=2)
    assert sim.cols[0] == "02_1"
    assert sim.cols[1] == "02_2"
    assert sim.cols[2] == "04_1"
    assert len(sim.cols) == 24


def test_reproducible_with_same_seed():
    sim1 = simulate(nrows=50, rseed=123)
    sim2 = simulate(nrows=50, rseed=123)
    assert np.array_equal(sim1.circ, sim2.circ)


def test_different_seeds_differ():
    sim1 = simulate(nrows=50, rseed=123)
    sim2 = simulate(nrows=50, rseed=456)
    assert not np.array_equal(sim1.circ, sim2.circ)


def test_pcirc_produces_expected_fraction():
    sim = simulate(nrows=1000, pcirc=0.5, rseed=42)
    frac = sim.circ.mean()
    assert 0.35 < frac < 0.65


def test_generate_pool_map_creates_file(tmp_path):
    sim = simulate(tpoints=12, nrows=20, nreps=2)
    out = str(tmp_path / "pool_map")
    sim.generate_pool_map(out_name=out)
    assert os.path.exists(out + ".parquet")
    pool_map = pd.read_parquet(out + ".parquet")["pool_number"].to_dict()
    assert isinstance(pool_map, dict)
    assert set(pool_map.keys()) == set(sim.cols)


def test_write_output_creates_files(tmp_path):
    sim = simulate(tpoints=12, nrows=20, nreps=2, rseed=42)
    out = str(tmp_path / "sim")
    sim.write_output(out_name=out)
    assert os.path.exists(out + "_with_noise.txt")
    assert os.path.exists(out + "_baseline.txt")
    assert os.path.exists(out + "_true_classes.txt")


def test_write_output_true_classes_shape(tmp_path):
    nrows = 30
    sim = simulate(tpoints=12, nrows=nrows, nreps=2, rseed=42)
    out = str(tmp_path / "sim")
    sim.write_output(out_name=out)
    df = pd.read_csv(out + "_true_classes.txt", sep="\t", index_col=0)
    assert len(df) == nrows
    assert "Circadian" in df.columns


def test_write_output_classes_binary(tmp_path):
    sim = simulate(tpoints=12, nrows=50, nreps=2, rseed=42)
    out = str(tmp_path / "sim")
    sim.write_output(out_name=out)
    df = pd.read_csv(out + "_true_classes.txt", sep="\t", index_col=0)
    assert set(df["Circadian"].unique()).issubset({0, 1})
