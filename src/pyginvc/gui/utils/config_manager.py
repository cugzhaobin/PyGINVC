"""
Config Manager - YAML configuration loading, saving, and validation.
"""
import os
import yaml


def load_yaml(path: str) -> dict:
    """Load a YAML configuration file and return the parsed dict."""
    with open(path, "r") as fid:
        cfg = yaml.safe_load(fid)
    return cfg


def save_yaml(cfg: dict, path: str):
    """Save a configuration dict to a YAML file."""
    with open(path, "w") as fid:
        yaml.dump(cfg, fid, default_flow_style=False, sort_keys=False)


def cfg_to_yaml_string(cfg: dict) -> str:
    """Convert a configuration dict to a YAML string."""
    return yaml.dump(cfg, default_flow_style=False, sort_keys=False)


def validate_keys(cfg: dict) -> list:
    """Check that all required top-level keys exist in the config."""
    required = [
        "dict_data",
        "dict_fault",
        "dict_green",
        "dict_weight",
        "dict_bound",
        "dict_export",
    ]
    missing = [k for k in required if k not in cfg]
    return missing
