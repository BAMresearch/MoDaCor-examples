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

    assert len(spec["nodes"]) == 192
    assert sum("origin" in node for node in spec["nodes"]) == 180
    assert Counter(node["module"] for node in spec["nodes"]) == {
        "AppendProcessingData": 36,
        "CopyDataBundleKeys": 27,
        "ThresholdMask": 21,
        "ApplyMask": 13,
        "Divide": 12,
        "IndexedAverager": 11,
        "DivideDatabundles": 18,
        "ConcatenateDatabundles": 10,
        "MultiplyDatabundles": 10,
        "AngleToQ": 8,
        "FindScaleFactor1D": 6,
        "FindCenterOfMass1D": 4,
        "IndexByCoordinate": 3,
        "Integrate1D": 2,
        "SubtractInterpolated1D": 1,
        "BitwiseOrMasks": 9,
        "Negate": 1,
    }
    assert "step_blocks" in pipeline.authored_spec
    assert "step_blocks" not in yaml.safe_load(pipeline.to_yaml())

    node_configs = {node["id"]: node["config"] for node in spec["nodes"]}
    assert node_configs["IP"]["bin_min"] == 0.002
    expected_noise_multipliers = {
        "SLF": 5,
        "SLR": 7,
        "SHF": 0,
        "SHR": 7,
        "BLF": 5,
        "BLR": 7,
        "BHF": 0,
        "BHR": 7,
    }
    for readout, multiplier in expected_noise_multipliers.items():
        copy_config = node_configs[f"prepare_diode.{readout}.copy_signal_to_dark_noise"]
        divide_config = node_configs[f"prepare_diode.{readout}.calculate_signal_to_dark_noise"]
        noise_config = node_configs[f"prepare_diode.{readout}.noise_mask"]
        assert copy_config["with_processing_keys"] == [readout, readout]
        assert copy_config["key_map"] == {"signal": "signal_to_dark_noise"}
        assert divide_config["with_processing_keys"] == [readout, readout]
        assert divide_config["dividend_data_key"] == "signal_to_dark_noise"
        assert divide_config["divisor_data_key"] == "dark_noise_std"
        assert noise_config["source_basedata_key"] == "signal_to_dark_noise"
        assert noise_config["lower_bound"] == multiplier
    for acquisition in ("S", "B"):
        for gain in ("low", "high"):
            center_pool = node_configs[f"prepare_acquisition_center.{acquisition}.pool_{gain}"]
            center_average = node_configs[f"prepare_acquisition_center.{acquisition}.average_{gain}"]
            science_pool = node_configs[f"combine_acquisition_diode_pairs.{acquisition}.pool_{gain}"]
            science_average = node_configs[f"combine_acquisition_diode_pairs.{acquisition}.average_{gain}"]
            assert center_pool["source_position_key"] == "pair_index"
            assert center_pool["alignment_key"] == "yaw"
            assert science_pool["source_position_key"] == "pair_index"
            assert science_pool["alignment_key"] == "Q"
            assert center_average["index_key"] == "pair_index"
            assert science_average["index_key"] == "pair_index"
            assert center_average["uncertainty_weight_key"] == "subread_sem"
            assert science_average["uncertainty_weight_key"] == "subread_sem"

        diode_fit = node_configs[f"combine_acquisition_diode_pairs.{acquisition}.fit_low_rear_to_front"]
        gain_fit = node_configs[f"scale_gain.{acquisition}.fit"]
        assert (diode_fit["fit_min_val"], diode_fit["fit_max_val"]) == (-0.0008, 0.0008)
        assert (gain_fit["fit_min_val"], gain_fit["fit_max_val"]) == (-0.005, -0.002)
        assert diode_fit["scale_uncertainty_key"] == "diode_scale_fit"
        assert gain_fit["scale_uncertainty_key"] == "gain_scale_fit"

    assert node_configs["scale_gain.S.fit"]["with_processing_keys"] == ["SH", "SL"]
    assert node_configs["scale_gain.B.fit"]["with_processing_keys"] == ["BH", "BL"]
    for item in ("sample", "background"):
        pool_config = node_configs[f"merge_acquisition.{item}.pool"]
        assert pool_config["uncertainty_key_policy"] == "fill_zero"


def test_usaxs_pipeline_graphs_group_each_authored_step_block() -> None:
    pipeline = Pipeline.from_yaml_file(PIPELINE_PATH)

    dot_source = pipeline.to_dot()
    mermaid_source = pipeline.to_mermaid(direction="TD")
    for block_id in pipeline.authored_spec["step_blocks"]:
        assert f'label="{block_id} (for_each)"' in dot_source
        assert f'["{block_id} (for_each)"]' in mermaid_source

    assert "subgraph" not in pipeline.to_dot(group_step_blocks=False)
    assert "subgraph" not in pipeline.to_mermaid(group_step_blocks=False)
