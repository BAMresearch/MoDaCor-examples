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

    assert len(spec["nodes"]) == 155
    assert sum("origin" in node for node in spec["nodes"]) == 152
    assert Counter(node["module"] for node in spec["nodes"]) == {
        "AppendProcessingData": 28,
        "ThresholdMask": 16,
        "ApplyMask": 16,
        "CopyDataBundleKeys": 16,
        "DivideDatabundles": 14,
        "Divide": 12,
        "AngleToQ": 8,
        "FindScaleFactor1D": 7,
        "MultiplyDatabundles": 7,
        "ConcatenateDatabundles": 5,
        "IndexedAverager": 5,
        "FindCenterOfMass1D": 4,
        "Integrate1D": 4,
        "SubtractInterpolated1D": 4,
        "BitwiseOrMasks": 4,
        "Negate": 4,
        "IndexByCoordinate": 1,
    }
    assert "step_blocks" in pipeline.authored_spec
    assert "step_blocks" not in yaml.safe_load(pipeline.to_yaml())

    node_configs = {node["id"]: node["config"] for node in spec["nodes"]}
    assert node_configs["IP"]["bin_min"] == 0.002
    assert {
        node_configs["scale_readout.SLR.fit"]["fit_min_val"],
        node_configs["scale_readout.SHF.fit"]["fit_min_val"],
        node_configs["scale_readout.SHR.fit"]["fit_min_val"],
    } == {0.002}
    for pair in ("SL", "SH", "BL", "BH"):
        pool_config = node_configs[f"prepare_pair_center.{pair}.pool"]
        average_config = node_configs[f"prepare_pair_center.{pair}.average"]
        assert pool_config["source_position_key"] == "pair_index"
        assert pool_config["alignment_key"] == "yaw"
        assert average_config["index_key"] == "pair_index"
        assert average_config["uncertainty_weight_key"] == "subread_sem"


def test_usaxs_pipeline_graphs_group_each_authored_step_block() -> None:
    pipeline = Pipeline.from_yaml_file(PIPELINE_PATH)

    dot_source = pipeline.to_dot()
    mermaid_source = pipeline.to_mermaid(direction="TD")
    for block_id in pipeline.authored_spec["step_blocks"]:
        assert f'label="{block_id} (for_each)"' in dot_source
        assert f'["{block_id} (for_each)"]' in mermaid_source

    assert "subgraph" not in pipeline.to_dot(group_step_blocks=False)
    assert "subgraph" not in pipeline.to_mermaid(group_step_blocks=False)
