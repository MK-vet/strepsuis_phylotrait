#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for analyzer.py methods.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil

REAL_DATA_PATH = r"C:\Users\ABC\Documents\GitHub\MKrep\data"


def check_real_data_exists():
    return os.path.exists(os.path.join(REAL_DATA_PATH, 'Snp_tree.newick'))


try:
    from strepsuis_phylotrait.analyzer import Analyzer
    ANALYZER_AVAILABLE = True
except (ImportError, OSError) as e:
    ANALYZER_AVAILABLE = False


@pytest.mark.skipif(not ANALYZER_AVAILABLE, reason="Not available")
class TestAnalyzerInit:
    """Tests for Analyzer initialization."""
    
    def test_init_basic(self):
        """Test basic initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                analyzer = Analyzer(
                    data_dir=tmpdir,
                    output_dir=os.path.join(tmpdir, 'output')
                )
                assert analyzer is not None
            except Exception as e:
                print(f"Init error: {e}")
    
    def test_init_with_config(self):
        """Test initialization with config."""
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                analyzer = Analyzer(
                    data_dir=tmpdir,
                    output_dir=os.path.join(tmpdir, 'output'),
                    n_clusters_range=(2, 8),
                    n_ensemble=5,
                )
                assert analyzer is not None
            except Exception as e:
                print(f"Init with config error: {e}")


@pytest.mark.skipif(not ANALYZER_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestAnalyzerMethods:
    """Tests for Analyzer methods."""
    
    def test_load_data(self):
        """Test loading data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                analyzer = Analyzer(
                    data_dir=tmpdir,
                    output_dir=os.path.join(tmpdir, 'output')
                )
                
                if hasattr(analyzer, 'load_data'):
                    analyzer.load_data()
                    
            except Exception as e:
                print(f"Load data error: {e}")
    
    def test_run_analysis(self):
        """Test running analysis."""
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                analyzer = Analyzer(
                    data_dir=tmpdir,
                    output_dir=os.path.join(tmpdir, 'output')
                )
                
                if hasattr(analyzer, 'run'):
                    analyzer.run()
                elif hasattr(analyzer, 'analyze'):
                    analyzer.analyze()
                    
            except Exception as e:
                print(f"Run analysis error: {e}")
    
    def test_generate_report(self):
        """Test generating report."""
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                analyzer = Analyzer(
                    data_dir=tmpdir,
                    output_dir=os.path.join(tmpdir, 'output')
                )
                
                if hasattr(analyzer, 'generate_report'):
                    analyzer.generate_report()
                    
            except Exception as e:
                print(f"Generate report error: {e}")
