#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for Visualizer and TraitAnalyzer classes.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np

try:
    from Bio import Phylo
    from io import StringIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        Visualizer,
        TraitAnalyzer,
        MCAAnalyzer,
        HTMLReportGenerator,
    )
    CLASSES_AVAILABLE = True
except (ImportError, OSError):
    CLASSES_AVAILABLE = False


class TestVisualizer:
    """Test Visualizer class."""
    
    @pytest.fixture
    def temp_output_folder(self):
        """Create temporary output folder."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    @pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
    def test_visualizer_init(self, temp_output_folder):
        """Test Visualizer initialization."""
        viz = Visualizer(temp_output_folder)
        
        assert viz is not None
        assert viz.output_folder == temp_output_folder


class TestTraitAnalyzer:
    """Test TraitAnalyzer class."""
    
    @pytest.fixture
    def sample_trait_data(self):
        """Create sample trait data."""
        return pd.DataFrame({
            'Trait_A': [1, 0, 1, 0, 1, 0],
            'Trait_B': [0, 1, 0, 1, 0, 1],
            'Trait_C': [1, 1, 0, 0, 1, 1],
        }, index=['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
    
    @pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
    def test_trait_analyzer_init(self, sample_trait_data):
        """Test TraitAnalyzer initialization."""
        analyzer = TraitAnalyzer(sample_trait_data)
        
        assert analyzer is not None


class TestMCAAnalyzer:
    """Test MCAAnalyzer class."""
    
    @pytest.fixture
    def sample_data(self):
        """Create sample data for MCA."""
        return pd.DataFrame({
            'A': [1, 0, 1, 0, 1, 0],
            'B': [0, 1, 0, 1, 0, 1],
            'C': [1, 1, 0, 0, 1, 1],
        })
    
    @pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
    def test_mca_analyzer_init(self, sample_data):
        """Test MCAAnalyzer initialization."""
        mca = MCAAnalyzer(sample_data)
        
        assert mca is not None


class TestHTMLReportGenerator:
    """Test HTMLReportGenerator class."""
    
    @pytest.fixture
    def temp_output_folder(self):
        """Create temporary output folder."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    @pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
    def test_html_report_generator_init(self, temp_output_folder):
        """Test HTMLReportGenerator initialization."""
        generator = HTMLReportGenerator(temp_output_folder)
        
        assert generator is not None
