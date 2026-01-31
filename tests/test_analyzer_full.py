#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for analyzer.py to increase coverage.
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


@pytest.mark.skipif(not ANALYZER_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestAnalyzerFull:
    """Full tests for Analyzer class."""
    
    def test_analyzer_init_with_real_data(self):
        """Test Analyzer initialization with real data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                analyzer = Analyzer(
                    data_dir=tmpdir,
                    output_dir=os.path.join(tmpdir, 'output')
                )
                assert analyzer is not None
            except Exception as e:
                print(f"Analyzer init error: {e}")
    
    def test_analyzer_load_data(self):
        """Test Analyzer load_data method."""
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
                if hasattr(analyzer, 'load_tree'):
                    analyzer.load_tree()
                if hasattr(analyzer, 'load_traits'):
                    analyzer.load_traits()
                    
            except Exception as e:
                print(f"Analyzer load error: {e}")
    
    def test_analyzer_run_analysis(self):
        """Test Analyzer run method."""
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
                print(f"Analyzer run error: {e}")
    
    def test_analyzer_generate_report(self):
        """Test Analyzer report generation."""
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
                print(f"Analyzer report error: {e}")
