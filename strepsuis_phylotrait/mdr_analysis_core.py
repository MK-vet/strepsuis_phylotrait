"""Compatibility stub for tests that patch mdr_analysis_core.setup_environment."""

from typing import Optional


def setup_environment(csv_path: Optional[str] = None, output_folder: str = "output") -> Optional[str]:
    """Return the provided csv_path for test compatibility."""
    return csv_path
