"""
Tests for strepsuis_phylotrait

Basic test suite to ensure the package is functional.
"""

import pytest

from strepsuis_phylotrait import __version__


def test_version():
    """Test that version is defined."""
    assert __version__ == "1.0.0"


def test_imports():
    """Test that main classes can be imported."""
    from strepsuis_phylotrait import Config, PhyloTraitAnalyzer

    assert PhyloTraitAnalyzer is not None
    assert Config is not None


def test_config_initialization():
    """Test Config class initialization."""
    import tempfile
    from strepsuis_phylotrait import Config

    with tempfile.TemporaryDirectory() as tmpdir:
        config = Config(data_dir=tmpdir, output_dir=tmpdir)

        assert config.data_dir == tmpdir
        assert config.output_dir == tmpdir


def test_analyzer_initialization():
    """Test PhyloTraitAnalyzer initialization."""
    import tempfile
    from strepsuis_phylotrait import PhyloTraitAnalyzer

    with tempfile.TemporaryDirectory() as tmpdir:
        analyzer = PhyloTraitAnalyzer(data_dir=tmpdir, output_dir=tmpdir)

        assert analyzer.data_dir == tmpdir
        assert analyzer.output_dir == tmpdir


@pytest.mark.integration
def test_example_data_exists():
    """Test that example data files exist."""
    from pathlib import Path

    examples_dir = Path("examples")
    if examples_dir.exists():
        csv_files = list(examples_dir.glob("*.csv"))
        assert len(csv_files) > 0, "No CSV files found in examples directory"
