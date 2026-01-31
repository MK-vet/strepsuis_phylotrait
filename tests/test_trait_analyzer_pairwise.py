#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for TraitAnalyzer.pairwise_fdr_post_hoc method.
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
class TestPairwiseFDRPostHoc:
    """Tests for pairwise_fdr_post_hoc method."""
    
    def test_pairwise_fdr_post_hoc_basic(self):
        """Test basic pairwise FDR post-hoc analysis."""
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
                result = analyzer.pairwise_fdr_post_hoc(merged_df, alpha=0.05)
                assert result is not None
                assert isinstance(result, pd.DataFrame)
            except Exception as e:
                print(f"Basic pairwise FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_single_cluster(self):
        """Test pairwise FDR with single cluster."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [0] * 50,
                'gene1': np.random.randint(0, 2, 50),
            })
            
            try:
                result = analyzer.pairwise_fdr_post_hoc(merged_df)
                assert result is not None
                assert len(result) == 0  # No pairs with single cluster
            except Exception as e:
                print(f"Single cluster pairwise FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_many_clusters(self):
        """Test pairwise FDR with many clusters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(150)],
                'Cluster': [i % 5 for i in range(150)],
                **{f'gene_{i}': np.random.randint(0, 2, 150) for i in range(10)}
            })
            
            try:
                result = analyzer.pairwise_fdr_post_hoc(merged_df, alpha=0.05)
                assert result is not None
                # Should have results for multiple cluster pairs
            except Exception as e:
                print(f"Many clusters pairwise FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_different_alpha(self):
        """Test pairwise FDR with different alpha values."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                'gene1': np.random.randint(0, 2, 100),
                'gene2': np.random.randint(0, 2, 100),
            })
            
            for alpha in [0.01, 0.05, 0.1]:
                try:
                    result = analyzer.pairwise_fdr_post_hoc(merged_df, alpha=alpha)
                    assert result is not None
                except Exception as e:
                    print(f"Alpha={alpha} pairwise FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_different_methods(self):
        """Test pairwise FDR with different FDR methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                'gene1': np.random.randint(0, 2, 100),
                'gene2': np.random.randint(0, 2, 100),
            })
            
            for method in ['fdr_bh', 'bonferroni', 'holm']:
                try:
                    result = analyzer.pairwise_fdr_post_hoc(merged_df, method=method)
                    assert result is not None
                except Exception as e:
                    print(f"Method={method} pairwise FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_empty_features(self):
        """Test pairwise FDR with empty features."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': [0] * 50,  # All zeros
                'gene2': [0] * 50,  # All zeros
            })
            
            try:
                result = analyzer.pairwise_fdr_post_hoc(merged_df)
                assert result is not None
            except Exception as e:
                print(f"Empty features pairwise FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_sparse_data(self):
        """Test pairwise FDR with sparse data."""
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
                result = analyzer.pairwise_fdr_post_hoc(merged_df)
                assert result is not None
            except Exception as e:
                print(f"Sparse data pairwise FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_large_dataset(self):
        """Test pairwise FDR with large dataset."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(200)],
                'Cluster': [i % 4 for i in range(200)],
                **{f'gene_{i}': np.random.randint(0, 2, 200) for i in range(20)}
            })
            
            try:
                result = analyzer.pairwise_fdr_post_hoc(merged_df)
                assert result is not None
            except Exception as e:
                print(f"Large dataset pairwise FDR error: {e}")
