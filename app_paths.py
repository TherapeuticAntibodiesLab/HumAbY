from __future__ import annotations

import os
from pathlib import Path

APP_SLUG = "antibody-humanizer"


def _xdg_dir(env_name: str, fallback: str) -> Path:
    value = os.environ.get(env_name)
    if value:
        return Path(value) / APP_SLUG
    return Path.home() / fallback / APP_SLUG


def get_data_dir() -> Path:
    return _xdg_dir("XDG_DATA_HOME", ".local/share")


def get_cache_dir() -> Path:
    return _xdg_dir("XDG_CACHE_HOME", ".cache")


def get_config_dir() -> Path:
    return _xdg_dir("XDG_CONFIG_HOME", ".config")