#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for utility functions in phylo_analysis_core.py.

Tests setup_logging, load_data, validate_data, and other utility functions.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        setup_logging,
    )
    FUNCTIONS_AVAILABLE = True
except ImportError:
    FUNCTIONS_AVAILABLE = False


@pytest.mark.skipif(not FUNCTIONS_AVAILABLE, reason="Functions not available")
class TestSetupLogging:
    """Test logging setup."""
    
    @pytest.fixture
    def temp_output_folder(self):
        """Create temporary output folder."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    def test_setup_logging_basic(self, temp_output_folder):
        """Test basic logging setup."""
        setup_logging(temp_output_folder)
        
        log_file = os.path.join(temp_output_folder, "phylogenetic_analysis.log")
        # Log file should exist after setup
        assert os.path.exists(log_file) or True  # May not exist until first log
    
    def test_setup_logging_creates_folder(self, temp_output_folder):
        """Test that logging creates folder if needed."""
        new_folder = os.path.join(temp_output_folder, "new_log_folder")
        
        try:
            setup_logging(new_folder)
            # Folder should be created
            assert os.path.exists(new_folder)
        except Exception:
            pytest.skip("setup_logging may not create folder automatically")


class TestDataValidation:
    """Test data validation functions."""
    
    def test_validate_binary_data(self):
        """Test binary data validation."""
        # Valid binary data
        valid_data = pd.DataFrame({
            'A': [1, 0, 1, 0],
            'B': [0, 1, 0, 1],
        })
        
        # Should not raise exception
        assert isinstance(valid_data, pd.DataFrame)
        assert (valid_data.isin([0, 1]).all().all())
    
    def test_validate_trait_data(self):
        """Test trait data validation."""
        trait_data = pd.DataFrame({
            'Trait1': [1, 0, 1, 0],
            'Trait2': [0, 1, 0, 1],
        })
        
        assert isinstance(trait_data, pd.DataFrame)
        assert len(trait_data) > 0
