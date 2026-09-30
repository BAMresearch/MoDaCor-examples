"""Lean preprocessing for the I22 four-scan USAXS acquisition layout."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from attrs import define, field

PREPROCESSING_VERSION = "2026-09-29-i22-usaxs-v1"
FRONT_DIODE_CHANNEL = 3
REAR_DIODE_CHANNEL = 3
I0_CHANNEL = 6


@define(frozen=True, slots=True)
class USAXSScanSet:
    """The low-dark, low-scan, high-dark, and high-scan acquisition numbers."""

    low_dark: int = field(converter=int)
    low_scan: int = field(converter=int)
    high_dark: int = field(converter=int)
    high_scan: int = field(converter=int)

    @classmethod
    def consecutive(cls, low_dark: int) -> "USAXSScanSet":
        """Construct the I22 acquisition pattern from its first scan number."""

        return cls(low_dark, low_dark + 1, low_dark + 2, low_dark + 3)

    @property
    def scan_numbers(self) -> tuple[int, int, int, int]:
        return self.low_dark, self.low_scan, self.high_dark, self.high_scan


@define(frozen=True, slots=True)
class ReducedReadout:
    signal: np.ndarray
    subread_sem: np.ndarray
    count_time: np.ndarray


@define(frozen=True, slots=True)
class DarkRate:
    mean: float
    sem: float
    exposure_count: int
    subreads_per_exposure: int


def _decode(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray) and value.shape == ():
        return _decode(value.item())
    return str(value)


def _scan_path(data_dir: Path, scan_number: int) -> Path:
    return data_dir / f"i22-{scan_number}.nxs"


def _broadcast_count_time(values: Any, shape: tuple[int, ...], *, source: str) -> np.ndarray:
    count_time = np.asarray(values, dtype=float).squeeze()
    if count_time.size == 1:
        count_time = np.full(shape, float(count_time.reshape(-1)[0]), dtype=float)
    else:
        try:
            count_time = np.broadcast_to(count_time, shape).astype(float, copy=True)
        except ValueError as exc:
            raise ValueError(f"Count-time shape at {source} cannot broadcast to {shape}.") from exc
    if np.any(~np.isfinite(count_time)) or np.any(count_time <= 0.0):
        raise ValueError(f"Count times at {source} must be finite and positive.")
    return count_time


def _reduce_subreads(values: np.ndarray, *, source: str) -> tuple[np.ndarray, np.ndarray, int]:
    values = np.asarray(values, dtype=float)
    if values.ndim < 2:
        raise ValueError(f"Expected a subread dimension at {source}, got shape {values.shape}.")
    valid = np.isfinite(values)
    counts = np.sum(valid, axis=-1)
    if np.any(counts < 2):
        raise ValueError(f"Every exposure at {source} needs at least two finite subreads.")
    total = np.sum(np.where(valid, values, 0.0), axis=-1)
    mean = total / counts
    deviations = np.where(valid, values - mean[..., None], 0.0)
    variance = np.sum(deviations**2, axis=-1) / (counts - 1)
    sem = np.sqrt(variance / counts)
    return np.asarray(mean).squeeze(), np.asarray(sem).squeeze(), int(values.shape[-1])


def _read_scan_readout(source: h5py.File, data_path: str, count_time_path: str, channel: int) -> ReducedReadout:
    values = np.asarray(source[data_path][()], dtype=float)
    if values.ndim != 3 or channel >= values.shape[-1]:
        raise ValueError(f"Unexpected scan readout shape {values.shape} at {data_path}.")
    mean, sem, _ = _reduce_subreads(values[..., channel], source=data_path)
    if mean.ndim != 1:
        raise ValueError(f"Reduced scan readout at {data_path} is not one-dimensional: {mean.shape}.")
    count_time = _broadcast_count_time(source[count_time_path][()], mean.shape, source=count_time_path)
    return ReducedReadout(mean, sem, count_time)


def _read_dark_rate(source: h5py.File, data_path: str, count_time_path: str, channel: int) -> DarkRate:
    values = np.asarray(source[data_path][()], dtype=float)
    if values.ndim != 4 or channel >= values.shape[-1]:
        raise ValueError(f"Unexpected dark readout shape {values.shape} at {data_path}.")
    exposure_means, _, subreads = _reduce_subreads(values[..., channel], source=data_path)
    exposure_means = np.asarray(exposure_means, dtype=float).reshape(-1)
    count_time = _broadcast_count_time(source[count_time_path][()], exposure_means.shape, source=count_time_path)
    rates = exposure_means / count_time
    rates = rates[np.isfinite(rates)]
    if rates.size < 2:
        raise ValueError(f"Dark readout at {data_path} needs at least two finite exposures.")
    return DarkRate(
        mean=float(np.mean(rates)),
        sem=float(np.std(rates, ddof=1) / np.sqrt(rates.size)),
        exposure_count=int(rates.size),
        subreads_per_exposure=subreads,
    )


def _validate_scan_titles(paths: dict[str, Path]) -> dict[str, str]:
    titles: dict[str, str] = {}
    expected = {
        "low_dark": ("low gain", "dc scan"),
        "low_scan": ("low gain", "usaxs scan"),
        "high_dark": ("high gain", "dc scan"),
        "high_scan": ("high gain", "usaxs scan"),
    }
    for role, path in paths.items():
        with h5py.File(path, "r") as source:
            title_path = "/entry1/title" if "/entry1/title" in source else "/entry/sample/title"
            title = _decode(source[title_path][()])
        titles[role] = title
        title_lower = title.lower()
        if not all(fragment in title_lower for fragment in expected[role]):
            raise ValueError(f"Scan assigned as {role!r} has unexpected title {title!r}.")
    return titles


def _signature(paths: dict[str, Path], scan_set: USAXSScanSet) -> str:
    inputs = {
        role: {"path": path.name, "size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns}
        for role, path in paths.items()
    }
    payload = {
        "version": PREPROCESSING_VERSION,
        "scan_numbers": scan_set.scan_numbers,
        "channels": {"front": FRONT_DIODE_CHANNEL, "rear": REAR_DIODE_CHANNEL, "i0": I0_CHANNEL},
        "inputs": inputs,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _write_dataset(group: h5py.Group, name: str, values: Any, units: str, **attrs: Any) -> h5py.Dataset:
    dataset = group.create_dataset(name, data=values)
    dataset.attrs["units"] = units
    for key, value in attrs.items():
        dataset.attrs[key] = value
    return dataset


def _write_curve(
    entry: h5py.Group,
    *,
    gain: str,
    diode: str,
    scan_number: int,
    dark_scan_number: int,
    scan: ReducedReadout,
    dark: DarkRate,
    yaw: np.ndarray,
    energy: float,
    i0: ReducedReadout,
) -> None:
    name = f"{gain}_gain_{diode}"
    group = entry.create_group(name)
    group.attrs.update(
        NX_class="NXdata",
        signal="signal",
        axes="yaw",
        scan_number=scan_number,
        dark_scan_number=dark_scan_number,
        diode=diode,
        amplifier_gain=gain,
    )
    adjusted = scan.signal - dark.mean * scan.count_time
    dark_offset_sem = np.abs(scan.count_time) * dark.sem
    _write_dataset(group, "signal", adjusted, "count", long_name="Exposure-adjusted dark-subtracted diode")
    _write_dataset(group, "subread_sem", scan.subread_sem, "count")
    _write_dataset(group, "dark_offset_sem", dark_offset_sem, "count")
    yaw_dataset = _write_dataset(group, "yaw", yaw, "microradian")
    yaw_dataset.attrs["indices"] = 0
    _write_dataset(group, "count_time", scan.count_time, "s")
    _write_dataset(group, "energy", np.asarray(energy), "keV")
    _write_dataset(group, "i0", i0.signal, "count")
    _write_dataset(group, "i0_subread_sem", i0.subread_sem, "count")
    _write_dataset(group, "i0_count_time", i0.count_time, "s")
    _write_dataset(group, "dark_rate", np.asarray(dark.mean), "count/s")
    _write_dataset(group, "dark_rate_sem", np.asarray(dark.sem), "count/s")
    group.attrs["dark_exposure_count"] = dark.exposure_count
    group.attrs["dark_subreads_per_exposure"] = dark.subreads_per_exposure


def preprocess_usaxs_acquisition(
    data_dir: str | Path,
    scan_set: USAXSScanSet,
    output_dir: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Reduce and dark-correct one four-scan I22 USAXS acquisition."""

    data_dir = Path(data_dir).resolve()
    output_dir = Path(output_dir).resolve()
    paths = {
        role: _scan_path(data_dir, scan_number)
        for role, scan_number in zip(
            ("low_dark", "low_scan", "high_dark", "high_scan"),
            scan_set.scan_numbers,
            strict=True,
        )
    }
    missing = [path for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing I22 USAXS scan files: " + ", ".join(str(path) for path in missing))
    titles = _validate_scan_titles(paths)
    signature = _signature(paths, scan_set)
    output_file = output_dir / f"i22-{scan_set.low_dark}-{scan_set.high_scan}_usaxs_modacor.nxs"
    if output_file.exists() and not overwrite:
        try:
            with h5py.File(output_file, "r") as existing:
                if _decode(existing.attrs.get("preprocessing_signature", "")) == signature:
                    return output_file
        except OSError:
            pass

    curves: dict[tuple[str, str], tuple[ReducedReadout, DarkRate]] = {}
    scan_metadata: dict[str, tuple[np.ndarray, float, ReducedReadout]] = {}
    detector_paths = {
        "front": ("/entry/HG_USAXS_DIODE/data", "/entry/instrument/HG_USAXS_DIODE/count_time", FRONT_DIODE_CHANNEL),
        "rear": ("/entry/BSDIODES/data", "/entry/instrument/BSDIODES/count_time", REAR_DIODE_CHANNEL),
    }
    dark_paths = {
        "front": ("/entry1/user_tetramm/data", "/entry1/instrument/user_tetramm/count_time", FRONT_DIODE_CHANNEL),
        "rear": ("/entry1/bsdiodes/data", "/entry1/instrument/bsdiodes/count_time", REAR_DIODE_CHANNEL),
    }
    for gain in ("low", "high"):
        with h5py.File(paths[f"{gain}_dark"], "r") as dark_source, h5py.File(
            paths[f"{gain}_scan"], "r"
        ) as scan_source:
            for diode in ("front", "rear"):
                scan_readout = _read_scan_readout(scan_source, *detector_paths[diode])
                dark_rate = _read_dark_rate(dark_source, *dark_paths[diode])
                curves[(gain, diode)] = scan_readout, dark_rate
            i0 = _read_scan_readout(
                scan_source,
                "/entry/I0/data",
                "/entry/instrument/I0/count_time",
                I0_CHANNEL,
            )
            yaw = np.asarray(scan_source["/entry/instrument/downstream_yaw/value"][()], dtype=float).squeeze()
            energy = float(np.asarray(scan_source["/entry/instrument/monochromator/energy"][()]).squeeze())
            if yaw.ndim != 1 or yaw.shape != i0.signal.shape:
                raise ValueError(f"{gain}-gain yaw and I0 shapes differ: {yaw.shape} and {i0.signal.shape}.")
            for diode in ("front", "rear"):
                if curves[(gain, diode)][0].signal.shape != yaw.shape:
                    raise ValueError(f"{gain}-gain {diode} signal and yaw shapes differ.")
            scan_metadata[gain] = yaw, energy, i0

    output_dir.mkdir(parents=True, exist_ok=True)
    temporary = output_file.with_suffix(output_file.suffix + ".tmp")
    with h5py.File(temporary, "w") as target:
        target.attrs.update(
            creator="MoDaCor-examples I22 USAXS preprocessor",
            preprocessing_version=PREPROCESSING_VERSION,
            preprocessing_signature=signature,
        )
        entry = target.create_group("entry1")
        entry.attrs.update(NX_class="NXentry", default="low_gain_rear")
        entry.create_dataset("title", data=titles["low_scan"])
        provenance = entry.create_group("provenance")
        provenance.attrs["NX_class"] = "NXcollection"
        for role, scan_number in zip(paths, scan_set.scan_numbers, strict=True):
            provenance.create_dataset(f"{role}_scan_number", data=scan_number)
            provenance.create_dataset(f"{role}_source_file", data=paths[role].name)
            provenance.create_dataset(f"{role}_title", data=titles[role])
        provenance.create_dataset(
            "dark_correction",
            data="adjusted = scan_mean - mean(dark_exposure_mean / dark_exposure_time) * scan_time",
        )
        provenance.create_dataset("i0_dark_correction", data="not applied")

        for gain, scan_number, dark_scan_number in (
            ("low", scan_set.low_scan, scan_set.low_dark),
            ("high", scan_set.high_scan, scan_set.high_dark),
        ):
            yaw, energy, i0 = scan_metadata[gain]
            for diode in ("front", "rear"):
                scan, dark = curves[(gain, diode)]
                _write_curve(
                    entry,
                    gain=gain,
                    diode=diode,
                    scan_number=scan_number,
                    dark_scan_number=dark_scan_number,
                    scan=scan,
                    dark=dark,
                    yaw=yaw,
                    energy=energy,
                    i0=i0,
                )
    temporary.replace(output_file)
    return output_file


BACKGROUND_SCAN_SET = USAXSScanSet.consecutive(977724)
SAMPLE_SCAN_SETS = (
    USAXSScanSet.consecutive(978497),
    USAXSScanSet.consecutive(978502),
    USAXSScanSet.consecutive(978507),
)
