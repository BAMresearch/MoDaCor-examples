from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import yaml
from modacor.runner.pipeline import Pipeline

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from i22_helpers import (
    COMBINED_PIPELINES,
    WAXS_PIPELINES,
    preprocess_measurement,
    sample_aligned_paths,
    sample_static_paths,
)


def _write_measurement(
    path: Path,
    *,
    bsdiodes: np.ndarray,
    i0: np.ndarray,
) -> None:
    leading_shape = bsdiodes.shape[:2]
    with h5py.File(path, "w") as h5:
        entry = h5.create_group("entry1")
        entry.attrs["NX_class"] = "NXentry"
        entry.create_dataset("detector/data", data=np.zeros((*leading_shape, 2, 3)))
        entry.create_dataset("Pilatus2M_WAXS/data", data=np.zeros((*leading_shape, 3, 2)))
        entry.create_dataset("bsdiodes/data", data=bsdiodes)
        entry.create_dataset("I0/data", data=i0)
        entry.create_dataset("instrument/detector/count_time", data=np.ones(leading_shape))
        entry.create_dataset("instrument/Pilatus2M_WAXS/count_time", data=np.ones(leading_shape))
        sample = entry.create_group("sample")
        sample.attrs["NX_class"] = "NXsample"
        sample.create_dataset("thickness", data=0.5)


def _readout(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    channel = values[..., 1]
    return np.mean(channel, axis=-1), np.std(channel, axis=-1, ddof=1) / np.sqrt(channel.shape[-1])


def test_preprocessing_calibrates_transmission_and_propagates_sem(tmp_path: Path) -> None:
    measurement_bsdiodes = np.zeros((1, 2, 4, 2))
    measurement_i0 = np.zeros((1, 2, 4, 2))
    measurement_bsdiodes[..., 1] = [[(3.0, 4.0, 5.0, 4.0), (7.0, 8.0, 9.0, 8.0)]]
    measurement_i0[..., 1] = [[(1.0, 2.0, 3.0, 2.0), (1.0, 2.0, 3.0, 2.0)]]

    # Deliberately use a different leading shape to confirm that the reference
    # supplies one portable scalar calibration rather than frame-paired values.
    reference_bsdiodes = np.zeros((1, 3, 2, 2))
    reference_i0 = np.zeros((1, 3, 2, 2))
    reference_bsdiodes[..., 1] = [[[5.0, 7.0], [6.0, 6.0], [7.0, 5.0]]]
    reference_i0[..., 1] = [[[1.0, 3.0], [2.0, 2.0], [3.0, 1.0]]]

    measurement = tmp_path / "measurement.nxs"
    reference = tmp_path / "open_beam.nxs"
    _write_measurement(measurement, bsdiodes=measurement_bsdiodes, i0=measurement_i0)
    _write_measurement(reference, bsdiodes=reference_bsdiodes, i0=reference_i0)

    output = preprocess_measurement(
        measurement,
        tmp_path / "preprocessed",
        transmission_reference_file=reference,
    )

    diode_mean, diode_sem = _readout(measurement_bsdiodes)
    i0_mean, i0_sem = _readout(measurement_i0)
    reference_diode = reference_bsdiodes[..., 1].reshape(-1)
    reference_i0_values = reference_i0[..., 1].reshape(-1)
    reference_diode_mean = np.mean(reference_diode)
    reference_i0_mean = np.mean(reference_i0_values)
    reference_diode_sem = np.std(reference_diode, ddof=1) / np.sqrt(reference_diode.size)
    reference_i0_sem = np.std(reference_i0_values, ddof=1) / np.sqrt(reference_i0_values.size)
    reference_ratio = reference_diode_mean / reference_i0_mean
    reference_ratio_sem = reference_ratio * np.hypot(
        reference_diode_sem / reference_diode_mean,
        reference_i0_sem / reference_i0_mean,
    )
    expected = (diode_mean / i0_mean) / reference_ratio
    expected_sem = expected * np.sqrt(
        (diode_sem / diode_mean) ** 2
        + (i0_sem / i0_mean) ** 2
        + (reference_ratio_sem / reference_ratio) ** 2
    )

    with h5py.File(output, "r") as h5:
        np.testing.assert_allclose(h5["/entry1/sample/transmission"][..., 0, 0], expected)
        np.testing.assert_allclose(h5["/entry1/sample/transmission_sem"][..., 0, 0], expected_sem)
        np.testing.assert_allclose(h5["/modacor/normalization/i0_channel_1_mean"][..., 0, 0], i0_mean)
        np.testing.assert_allclose(h5["/modacor/calibration/bsdiodes_to_i0_ratio"][()], reference_ratio)
        np.testing.assert_allclose(h5["/modacor/calibration/bsdiodes_to_i0_ratio_sem"][()], reference_ratio_sem)
        assert h5["/entry1/sample/transmission"].attrs["units"] == "dimensionless"
        assert h5["/entry1/sample/transmission_sem"].attrs["units"] == "dimensionless"
        assert h5["/modacor/normalization/i0_channel_1_mean"].attrs["units"] == "count/s"
        assert h5["/modacor/normalization/i0_channel_1_sem"].attrs["units"] == "count/s"
        assert h5["/modacor/normalization/bsdiodes_channel_1_mean"].attrs["units"] == "count/s"
        assert h5["/modacor/calibration/reference_i0_mean"].attrs["units"] == "count/s"
        assert h5["/modacor/calibration/reference_bsdiodes_mean"].attrs["units"] == "count/s"
        assert h5["/modacor/normalization"].attrs["readout_units"] == "count/s"
        assert isinstance(h5.get("/entry1/sample/thickness", getlink=True), h5py.ExternalLink)

    with h5py.File(measurement, "r") as source:
        assert "/entry1/sample/transmission" not in source


def test_reference_file_changes_invalidate_cached_preprocessing(tmp_path: Path) -> None:
    values = np.zeros((1, 1, 3, 2))
    values[..., 1] = 2.0
    measurement = tmp_path / "measurement.nxs"
    reference_one = tmp_path / "open_beam_one.nxs"
    reference_two = tmp_path / "open_beam_two.nxs"
    _write_measurement(measurement, bsdiodes=values * 2.0, i0=values)
    _write_measurement(reference_one, bsdiodes=values * 2.0, i0=values)
    _write_measurement(reference_two, bsdiodes=values * 4.0, i0=values)

    output = preprocess_measurement(
        measurement,
        tmp_path / "preprocessed",
        transmission_reference_file=reference_one,
    )
    with h5py.File(output, "r") as h5:
        np.testing.assert_allclose(h5["/entry1/sample/transmission"][()], 1.0)

    preprocess_measurement(
        measurement,
        tmp_path / "preprocessed",
        transmission_reference_file=reference_two,
    )
    with h5py.File(output, "r") as h5:
        np.testing.assert_allclose(h5["/entry1/sample/transmission"][()], 0.5)
        assert h5["/modacor/calibration"].attrs["transmission_reference_file"].endswith(
            "open_beam_two.nxs"
        )


def test_preprocessing_can_record_known_calibration_sample_thickness(tmp_path: Path) -> None:
    values = np.zeros((1, 1, 3, 2))
    values[..., 1] = 2.0
    measurement = tmp_path / "measurement.nxs"
    reference = tmp_path / "open_beam.nxs"
    _write_measurement(measurement, bsdiodes=values * 2.0, i0=values)
    _write_measurement(reference, bsdiodes=values * 2.0, i0=values)

    output = preprocess_measurement(
        measurement,
        tmp_path / "preprocessed",
        transmission_reference_file=reference,
        sample_thickness=1.0,
        sample_thickness_units="mm",
    )

    with h5py.File(output, "r") as h5:
        thickness = h5["/entry1/sample/thickness"]
        assert not isinstance(h5.get("/entry1/sample/thickness", getlink=True), h5py.ExternalLink)
        assert thickness[()] == 1.0
        assert thickness.attrs["units"] == "mm"
        assert thickness.attrs["source"] == "Known glassy-carbon calibration-sample thickness"


def test_operational_pipelines_normalize_time_before_transmission_and_flux() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    pipeline_files = {
        "SAXS": "I22_SAXS_solids_operando.yaml",
        **{f"WAXS_{profile}": filename for profile, filename in WAXS_PIPELINES.items()},
    }
    for detector_profile, filename in pipeline_files.items():
        detector = detector_profile.split("_", maxsplit=1)[0]
        pipeline = yaml.safe_load(
            (project_dir / "pipelines" / filename).read_text()
        )
        steps = pipeline["steps"]
        assert not any(step_id.startswith("BS_") for step_id in steps)

        for processing_key in ("sample", "background", "intensity_calibration"):
            suffix = f"_{processing_key}"
            assert steps[f"TI{suffix}"]["requires_steps"] == [f"PU{suffix}"]
            assert steps[f"TR{suffix}"]["requires_steps"] == [f"TI{suffix}"]
            assert steps[f"FL{suffix}"]["requires_steps"] == [f"TR{suffix}"]
            assert steps[f"FA{suffix}"]["requires_steps"] == [f"FL{suffix}"]

            source_ref = processing_key
            transmission = steps[f"TR{suffix}"]["configuration"]
            flux = steps[f"FL{suffix}"]["configuration"]
            assert transmission["divisor_source"] == f"{source_ref}::/entry1/sample/transmission"
            assert transmission["divisor_uncertainties_sources"] == {
                "transmission_SEM": f"{source_ref}::/entry1/sample/transmission_sem"
            }
            assert flux["divisor_source"] == (
                f"{source_ref}::/modacor/normalization/i0_channel_1_mean"
            )
            assert flux["divisor_uncertainties_sources"] == {
                "I0_SEM": f"{source_ref}::/modacor/normalization/i0_channel_1_sem"
            }

        assert steps["SP"]["configuration"]["with_processing_keys"] == ["sample"]
        assert steps["SP_intensity_calibration"]["configuration"]["with_processing_keys"] == [
            "intensity_calibration"
        ]
        assert steps["AV"]["configuration"]["with_processing_keys"] == ["sample"]
        assert steps["AV_intensity_calibration"]["configuration"]["with_processing_keys"] == [
            "intensity_calibration"
        ]
        assert steps["TH_intensity_calibration"]["configuration"]["divisor_source"] == (
            "intensity_calibration::/entry1/sample/thickness"
        )
        assert steps["TH_intensity_calibration"]["requires_steps"] == [
            "PO_intensity_calibration"
        ]
        assert steps["TH_sample"]["configuration"]["divisor_source"] == (
            "sample::/entry1/sample/thickness"
        )
        assert steps["TH_sample"]["requires_steps"] == ["PO"]
        assert steps["AV_intensity_calibration"]["requires_steps"] == [
            "TH_intensity_calibration"
        ]
        scale_fit = steps["SF"]["configuration"]
        assert scale_fit["with_processing_keys"] == [
            "intensity_calibration_fit",
            "intensity_calibration_reference",
        ]
        assert scale_fit["fit_model"] == "lognormal"
        assert scale_fit["uncertainty_weight_key"] == "combined"
        assert steps["AU"]["module"] == "MultiplyDatabundles"
        assert steps["AU"]["configuration"] == {
            "with_processing_keys": ["sample", "intensity_calibration_fit"],
            "multiplicand_data_key": "signal",
            "multiplier_data_key": "absolute_intensity_scale_factor",
        }
        assert steps["AU"]["requires_steps"] == ["CU", "SF"]
        assert "UL_intensity_calibration" not in steps
        assert "SF_units" not in steps

        aligned = sample_aligned_paths(detector)
        assert "/modacor/normalization/i0_channel_1_mean" in aligned
        assert "/modacor/normalization/i0_channel_1_sem" in aligned
        assert "/modacor/normalization/bsdiodes_channel_1_mean" not in aligned
        assert "/entry1/sample/thickness" in sample_static_paths()


def test_waxs_nosecone_profiles_select_aluminium_correction() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    pipelines = {
        profile: yaml.safe_load((project_dir / "pipelines" / filename).read_text())
        for profile, filename in WAXS_PIPELINES.items()
    }

    usaxs_steps = pipelines["usaxs_saxs_waxs"]["steps"]
    assert "AL" not in usaxs_steps
    assert usaxs_steps["PO"]["requires_steps"] == ["DE"]

    standard_steps = pipelines["standard_saxs_waxs"]["steps"]
    assert standard_steps["AL"]["module"] == "AttenuatorPlateCorrection"
    assert standard_steps["AL"]["configuration"]["thickness"] == 0.1
    assert standard_steps["AL"]["configuration"]["thickness_units"] == "cm"
    assert standard_steps["AL"]["configuration"]["apply_as"] == "divide"
    assert standard_steps["PO"]["requires_steps"] == ["AL"]


def test_combined_pipeline_expands_readable_detector_lanes() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    pipeline = Pipeline.from_yaml_file(
        project_dir / "pipelines" / COMBINED_PIPELINES["usaxs_saxs_waxs"]
    )
    spec = pipeline.to_spec()
    nodes = {node["id"]: node for node in spec["nodes"]}
    blocks = pipeline.authored_spec["step_blocks"]

    assert "USAXS" not in pipeline.name
    assert len(nodes) == 96
    assert sum("origin" in node for node in nodes.values()) == 87
    assert list(blocks) == [
        "normalize_frames",
        "detector_geometry",
        "correct_detector",
        "integrate_sample",
    ]
    assert list(blocks["normalize_frames"]["for_each"]) == [
        "saxs_sample",
        "waxs_sample",
        "saxs_background",
        "waxs_background",
        "saxs_glassy_carbon",
    ]

    processing_keys = {
        "saxs_sample": "sample",
        "waxs_sample": "waxs_sample",
        "saxs_background": "background",
        "waxs_background": "waxs_background",
        "saxs_glassy_carbon": "intensity_calibration",
    }
    local_chain = (
        "load",
        "raw_mask",
        "poisson",
        "normalize_time",
        "normalize_transmission",
        "normalize_flux",
        "average_frames",
        "reduce_raw_mask",
    )
    for item, processing_key in processing_keys.items():
        prefix = f"normalize_frames.{item}"
        for upstream, downstream in zip(local_chain, local_chain[1:]):
            assert nodes[f"{prefix}.{downstream}"]["requires_steps"] == [
                f"{prefix}.{upstream}"
            ]
        assert nodes[f"{prefix}.poisson"]["config"]["with_processing_keys"] == [
            processing_key
        ]

    for item, processing_key in {"saxs": "sample", "waxs": "waxs_sample"}.items():
        correction_prefix = f"correct_detector.{item}"
        integration_prefix = f"integrate_sample.{item}"
        assert nodes[f"{correction_prefix}.normalize_thickness"]["config"]["with_processing_keys"] == [
            processing_key
        ]
        assert nodes[f"{correction_prefix}.normalize_thickness"]["requires_steps"] == [
            f"{correction_prefix}.polarization"
        ]
        assert nodes[f"{integration_prefix}.azimuthal_average"]["requires_steps"] == [
            f"{integration_prefix}.plot_2d",
            f"{integration_prefix}.save_2d",
        ]
        assert nodes[f"{integration_prefix}.apply_absolute_scale"]["config"]["with_processing_keys"] == [
            processing_key,
            "intensity_calibration_fit",
        ]
        assert set(nodes[f"{integration_prefix}.apply_absolute_scale"]["requires_steps"]) == {
            f"{integration_prefix}.combine_uncertainties",
            "SF",
        }

    assert nodes["correct_detector.saxs_glassy_carbon.normalize_thickness"]["config"][
        "divisor_source"
    ] == "intensity_calibration::/entry1/sample/thickness"
    assert nodes["AV_intensity_calibration"]["requires_steps"] == [
        "correct_detector.saxs_glassy_carbon.normalize_thickness"
    ]
    assert nodes["CP_intensity_calibration_fit"]["requires_steps"] == [
        "CU_intensity_calibration"
    ]

    assert nodes["PD_intensity_calibration_reference_I"]["config"]["units_override"] == (
        "1/(cm*sr)"
    )
    assert nodes["SF"]["config"]["fit_max_val"] == 2.0
    assert nodes["SF"]["requires_steps"] == [
        "CP_intensity_calibration_fit",
        "PD_intensity_calibration_reference_Q",
    ]
    assert nodes["CAT_IQ"]["requires_steps"] == [
        "integrate_sample.saxs.apply_absolute_scale",
        "integrate_sample.waxs.apply_absolute_scale",
    ]
    concatenate = nodes["CAT_IQ"]["config"]
    assert concatenate["with_processing_keys"] == ["sample", "waxs_sample"]
    assert concatenate["data_keys"] == ["signal", "Q"]
    assert concatenate["output_processing_key"] == "saxs_waxs_1d"
    assert concatenate["sort_by"] == "Q"
    assert concatenate["source_index_key"] == "detector_index"
    assert concatenate["uncertainty_key_policy"] == "require_matching"
    assert nodes["PL_IQ"]["config"]["x_path"] == "/saxs_waxs_1d/Q/signal"
    assert nodes["PL_IQ"]["config"]["y_path"] == "/saxs_waxs_1d/signal/signal"
    assert nodes["PL_IQ"]["config"]["x_axis_type"] == "log"
    assert nodes["PL_IQ"]["config"]["y_axis_type"] == "log"
