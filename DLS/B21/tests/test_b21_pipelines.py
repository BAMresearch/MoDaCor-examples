from __future__ import annotations

from pathlib import Path

import yaml


PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipelines"


def _steps(filename: str) -> dict:
    return yaml.safe_load((PIPELINE_DIR / filename).read_text())["steps"]


def test_prefilter_normalizes_only_static_inputs_to_data_rank() -> None:
    steps = _steps("B21_frame_prefilter.yaml")

    signal_reduction = steps["RD_static_signal"]
    assert signal_reduction["module"] == "ReduceDimensionality"
    assert signal_reduction["configuration"]["axes"] == "non_data"
    assert signal_reduction["configuration"]["with_processing_keys"] == ["static"]

    mask_reduction = steps["RM_static_mask"]
    assert mask_reduction["module"] == "ReduceMask"
    assert mask_reduction["configuration"]["axes"] == "non_data"
    assert mask_reduction["configuration"]["with_processing_keys"] == ["static"]

    assert steps["PC_static"]["requires_steps"] == ["RD_static_signal"]
    assert "RM_static_mask" in steps["CP_static_sample"]["requires_steps"]
    assert not any(
        step["module"] in {"ReduceDimensionality", "ReduceMask"}
        and step["configuration"].get("with_processing_keys") != ["static"]
        for step in steps.values()
    )


def test_quality_pipeline_keeps_frame_dimensions_distinct() -> None:
    steps = _steps("B21_frame_quality.yaml")
    assert not any(
        step["module"] in {"ReduceDimensionality", "ReduceMask"}
        for step in steps.values()
    )
