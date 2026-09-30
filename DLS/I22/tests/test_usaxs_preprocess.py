from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from usaxs_preprocess import USAXSScanSet, preprocess_usaxs_acquisition


def _write_dark(path: Path, title: str, front_means: np.ndarray, rear_means: np.ndarray) -> None:
    shape = (1, front_means.size, 4, 11)
    front = np.zeros(shape)
    rear = np.zeros(shape)
    for index, value in enumerate(front_means):
        front[0, index, :, 3] = value
    for index, value in enumerate(rear_means):
        rear[0, index, :, 3] = value
    with h5py.File(path, "w") as h5:
        entry = h5.create_group("entry1")
        entry.create_dataset("title", data=title)
        entry.create_dataset("user_tetramm/data", data=front)
        entry.create_dataset("bsdiodes/data", data=rear)
        entry.create_dataset(
            "instrument/user_tetramm/count_time", data=np.full((1, front_means.size), 0.5)
        )
        entry.create_dataset(
            "instrument/bsdiodes/count_time", data=np.full((1, rear_means.size), 0.5)
        )


def _write_scan(path: Path, title: str, front_level: float, rear_level: float) -> None:
    shape = (5, 4, 11)
    front = np.zeros(shape)
    rear = np.zeros(shape)
    i0 = np.zeros(shape)
    offsets = np.array([-1.5, -0.5, 0.5, 1.5])
    front[:, :, 3] = front_level + offsets
    rear[:, :, 3] = rear_level + offsets
    i0[:, :, 6] = 100.0 + offsets
    with h5py.File(path, "w") as h5:
        entry = h5.create_group("entry")
        entry.create_dataset("sample/title", data=title)
        entry.create_dataset("HG_USAXS_DIODE/data", data=front)
        entry.create_dataset("BSDIODES/data", data=rear)
        entry.create_dataset("I0/data", data=i0)
        entry.create_dataset("instrument/HG_USAXS_DIODE/count_time", data=0.25)
        entry.create_dataset("instrument/BSDIODES/count_time", data=0.25)
        entry.create_dataset("instrument/I0/count_time", data=0.25)
        entry.create_dataset("instrument/downstream_yaw/value", data=np.linspace(-2.0, 2.0, 5))
        energy = entry.create_dataset("instrument/monochromator/energy", data=14.0)
        energy.attrs["units"] = "keV"


def test_preprocess_usaxs_acquisition_produces_four_dark_adjusted_curves(tmp_path: Path) -> None:
    scans = USAXSScanSet.consecutive(100)
    _write_dark(tmp_path / "i22-100.nxs", "sample - low gain dc scan", np.array([2, 4, 6]), np.array([4, 6, 8]))
    _write_scan(tmp_path / "i22-101.nxs", "sample - low gain usaxs scan", 20.0, 30.0)
    _write_dark(tmp_path / "i22-102.nxs", "sample - high gain dc scan", np.array([1, 2, 3]), np.array([3, 4, 5]))
    _write_scan(tmp_path / "i22-103.nxs", "sample - high gain usaxs scan", 10.0, 15.0)

    output = preprocess_usaxs_acquisition(tmp_path, scans, tmp_path / "processed")

    with h5py.File(output, "r") as h5:
        assert set(h5["/entry1"]) >= {
            "low_gain_front",
            "low_gain_rear",
            "high_gain_front",
            "high_gain_rear",
            "provenance",
        }
        # Low-front dark rates are [4, 8, 12] count/s, with mean 8.
        np.testing.assert_allclose(h5["/entry1/low_gain_front/signal"], 20.0 - 8.0 * 0.25)
        np.testing.assert_allclose(h5["/entry1/low_gain_front/dark_rate"], 8.0)
        expected_dark_rate_sem = np.std([4.0, 8.0, 12.0], ddof=1) / np.sqrt(3)
        np.testing.assert_allclose(
            h5["/entry1/low_gain_front/dark_offset_sem"],
            expected_dark_rate_sem * 0.25,
        )
        expected_subread_sem = np.std([-1.5, -0.5, 0.5, 1.5], ddof=1) / 2.0
        np.testing.assert_allclose(
            h5["/entry1/low_gain_front/subread_sem"], expected_subread_sem
        )
        np.testing.assert_allclose(h5["/entry1/low_gain_front/i0"], 100.0)
        assert h5["/entry1/low_gain_front/signal"].attrs["units"] == "count"
        assert h5["/entry1/low_gain_front/yaw"].attrs["units"] == "microradian"
        assert h5["/entry1/provenance/i0_dark_correction"].asstr()[()] == "not applied"


def test_preprocess_usaxs_acquisition_reuses_matching_cached_output(tmp_path: Path) -> None:
    scans = USAXSScanSet.consecutive(200)
    _write_dark(tmp_path / "i22-200.nxs", "sample - low gain dc scan", np.array([2, 2]), np.array([2, 2]))
    _write_scan(tmp_path / "i22-201.nxs", "sample - low gain usaxs scan", 20.0, 30.0)
    _write_dark(tmp_path / "i22-202.nxs", "sample - high gain dc scan", np.array([2, 2]), np.array([2, 2]))
    _write_scan(tmp_path / "i22-203.nxs", "sample - high gain usaxs scan", 20.0, 30.0)

    first = preprocess_usaxs_acquisition(tmp_path, scans, tmp_path / "processed")
    first_mtime = first.stat().st_mtime_ns
    second = preprocess_usaxs_acquisition(tmp_path, scans, tmp_path / "processed")

    assert first == second
    assert second.stat().st_mtime_ns == first_mtime
