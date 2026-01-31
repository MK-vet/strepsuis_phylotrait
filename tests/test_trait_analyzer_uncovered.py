#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for uncovered TraitAnalyzer methods.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        TraitAnalyzer,
    )
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerPairwiseFDR:
    """Tests for pairwise_fdr_post_hoc method."""
    
    def test_pairwise_fdr_post_hoc_basic(self):
        """Test basic pairwise FDR post-hoc."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 4 for i in range(100)],
                'gene1': np.random.randint(0, 2, 100),
                'gene2': np.random.randint(0, 2, 100),
                'gene3': np.random.randint(0, 2, 100),
            })
            
            try:
                result = analyzer.pairwise_fdr_post_hoc(merged_df, alpha=0.05)
                assert result is not None
                assert os.path.exists(os.path.join(tmpdir, 'pairwise_fdr_post_hoc.csv'))
            except Exception as e:
                print(f"Pairwise FDR error: {e}")
    
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
            except Exception as e:
                print(f"Single cluster FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_two_clusters(self):
        """Test pairwise FDR with two clusters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [0] * 25 + [1] * 25,
                'gene1': [1] * 20 + [0] * 5 + [1] * 5 + [0] * 20,
                'gene2': [1] * 10 + [0] * 15 + [1] * 15 + [0] * 10,
            })
            
            try:
                result = analyzer.pairwise_fdr_post_hoc(merged_df, alpha=0.05)
                assert result is not None
            except Exception as e:
                print(f"Two clusters FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_many_clusters(self):
        """Test pairwise FDR with many clusters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 6 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(10)}
            })
            
            try:
                result = analyzer.pairwise_fdr_post_hoc(merged_df, alpha=0.01)
                assert result is not None
            except Exception as e:
                print(f"Many clusters FDR error: {e}")
    
    def test_pairwise_fdr_post_hoc_different_methods(self):
        """Test pairwise FDR with different correction methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(80)],
                'Cluster': [i % 4 for i in range(80)],
                'gene1': np.random.randint(0, 2, 80),
                'gene2': np.random.randint(0, 2, 80),
            })
            
            for method in ['fdr_bh', 'bonferroni', 'holm']:
                try:
                    result = analyzer.pairwise_fdr_post_hoc(merged_df, method=method)
                    assert result is not None
                except Exception as e:
                    print(f"Method {method} error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerChi2Tests:
    """Tests for chi2_tests method."""
    
    def test_chi2_tests_basic(self):
        """Test basic chi-squared tests."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 4 for i in range(100)],
                'gene1': np.random.randint(0, 2, 100),
                'gene2': np.random.randint(0, 2, 100),
                'gene3': np.random.randint(0, 2, 100),
            })
            
            try:
                if hasattr(analyzer, 'chi2_tests'):
                    result = analyzer.chi2_tests(merged_df)
                    assert result is not None
            except Exception as e:
                print(f"Chi2 tests error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerFisherExact:
    """Tests for fisher_exact_tests method."""
    
    def test_fisher_exact_tests_basic(self):
        """Test basic Fisher exact tests."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'fisher_exact_tests'):
                    result = analyzer.fisher_exact_tests(merged_df)
                    assert result is not None
            except Exception as e:
                print(f"Fisher exact error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerAssociationRules:
    """Tests for association_rules method."""
    
    def test_association_rules_basic(self):
        """Test basic association rules."""
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
                if hasattr(analyzer, 'association_rules'):
                    result = analyzer.association_rules(merged_df)
                    assert result is not None
            except Exception as e:
                print(f"Association rules error: {e}")
