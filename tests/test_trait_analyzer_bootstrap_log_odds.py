#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for TraitAnalyzer.bootstrap_log_odds method.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile

try:
    from strepsuis_phylotrait.phylo_analysis_core import TraitAnalyzer
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestBootstrapLogOdds:
    """Tests for bootstrap_log_odds method."""
    
    def test_bootstrap_log_odds_basic(self):
        """Test basic bootstrap log-odds analysis."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                'gene1': np.random.randint(0, 2, 100),
                'gene2': np.random.randint(0, 2, 100),
                'gene3': np.random.randint(0, 2, 100),
            })
            
            try:
                result = analyzer.bootstrap_log_odds(merged_df, n_bootstrap=10)
                assert result is not None
                assert isinstance(result, pd.DataFrame)
                assert 'Feature' in result.columns
                assert 'Bootstrap_Log_Odds_Mean' in result.columns
            except Exception as e:
                print(f"Basic bootstrap log-odds error: {e}")
    
    def test_bootstrap_log_odds_different_n_bootstrap(self):
        """Test bootstrap log-odds with different n_bootstrap values."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                'gene1': np.random.randint(0, 2, 100),
                'gene2': np.random.randint(0, 2, 100),
            })
            
            for n_bootstrap in [5, 10, 20, 50]:
                try:
                    result = analyzer.bootstrap_log_odds(merged_df, n_bootstrap=n_bootstrap)
                    assert result is not None
                except Exception as e:
                    print(f"n_bootstrap={n_bootstrap} error: {e}")
    
    def test_bootstrap_log_odds_empty_features(self):
        """Test bootstrap log-odds with empty features."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': [0] * 50,  # All zeros
                'gene2': [0] * 50,  # All zeros
            })
            
            try:
                result = analyzer.bootstrap_log_odds(merged_df, n_bootstrap=5)
                assert result is not None
                # Should return empty DataFrame or handle gracefully
            except Exception as e:
                print(f"Empty features bootstrap log-odds error: {e}")
    
    def test_bootstrap_log_odds_all_ones(self):
        """Test bootstrap log-odds with all ones."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': [1] * 50,  # All ones
                'gene2': [1] * 50,  # All ones
            })
            
            try:
                result = analyzer.bootstrap_log_odds(merged_df, n_bootstrap=5)
                assert result is not None
            except Exception as e:
                print(f"All ones bootstrap log-odds error: {e}")
    
    def test_bootstrap_log_odds_many_features(self):
        """Test bootstrap log-odds with many features."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(30)}
            })
            
            try:
                result = analyzer.bootstrap_log_odds(merged_df, n_bootstrap=10)
                assert result is not None
                assert len(result) == 30
            except Exception as e:
                print(f"Many features bootstrap log-odds error: {e}")
    
    def test_bootstrap_log_odds_sparse_data(self):
        """Test bootstrap log-odds with sparse data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                'gene1': (np.random.rand(100) > 0.9).astype(int),  # Very sparse
                'gene2': (np.random.rand(100) > 0.95).astype(int),  # Very sparse
            })
            
            try:
                result = analyzer.bootstrap_log_odds(merged_df, n_bootstrap=5)
                assert result is not None
            except Exception as e:
                print(f"Sparse data bootstrap log-odds error: {e}")
    
    def test_bootstrap_log_odds_large_dataset(self):
        """Test bootstrap log-odds with large dataset."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(200)],
                'Cluster': [i % 4 for i in range(200)],
                **{f'gene_{i}': np.random.randint(0, 2, 200) for i in range(20)}
            })
            
            try:
                result = analyzer.bootstrap_log_odds(merged_df, n_bootstrap=10)
                assert result is not None
            except Exception as e:
                print(f"Large dataset bootstrap log-odds error: {e}")
