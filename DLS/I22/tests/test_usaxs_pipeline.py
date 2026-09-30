from __future__ import annotations

from collections import Counter
from pathlib import Path

import yaml

from modacor.runner.pipeline import Pipeline


I22_DIR = Path(__file__).resolve().parents[1]
PIPELINE_PATH = I22_DIR / "pipelines" / "I22_USAXS.yaml"


def test_usaxs_pipeline_expands_to_the_complete_execution_graph() -> None:
    pipeline = Pipeline.from_yaml_file(PIPELINE_PATH)
    spec = pipeline.to_spec()

    assert len(spec["nodes"]) == 123
    assert sum("origin" in node for node in spec["nodes"]) == 120
    assert Counter(node["module"] for node in spec["nodes"]) == {
        "AppendProcessingData": 28,
        "ThresholdMask": 16,
        "DivideDatabundles": 14,
        "Divide": 12,
        "ApplyMask": 8,
        "CopyDataBundleKeys": 8,
        "YawToQ": 8,
        "FindCenterOfMass1D": 4,
        "Integrate1D": 4,
        "SubtractInterpolated1D": 4,
        "BitwiseOrMasks": 4,
        "Negate": 4,
        "FindScaleFactor1D": 3,
        "MultiplyDatabundles": 3,
        "ConcatenateDatabundles": 1,
        "IndexPixels": 1,
        "IndexedAverager": 1,
    }
    assert "step_blocks" in pipeline.authored_spec
    assert "step_blocks" not in yaml.safe_load(pipeline.to_yaml())
