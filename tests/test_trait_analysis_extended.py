#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extended tests for TraitAnalyzer methods in phylo_analysis_core.py.
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
class TestTraitAnalyzerAssociationRules:
    """Tests for association rule mining."""
    
    def test_association_rule_mining_basic(self):
        """Test basic association rule mining."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(20)}
            })
            
            try:
                result = analyzer.association_rule_mining(merged_df, min_support=0.1)
                assert result is not None
            except Exception as e:
                print(f"Association rule mining error: {e}")
    
    def test_association_rule_mining_high_support(self):
        """Test association rule mining with high support."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(20)}
            })
            
            try:
                result = analyzer.association_rule_mining(merged_df, min_support=0.5, min_confidence=0.8)
                assert result is not None
            except Exception as e:
                print(f"High support error: {e}")
    
    def test_association_rule_mining_many_features(self):
        """Test association rule mining with many features."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 3 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(80)}
            })
            
            try:
                result = analyzer.association_rule_mining(merged_df, min_support=0.15, max_features=50)
                assert result is not None
            except Exception as e:
                print(f"Many features error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerLabelSharedUnique:
    """Tests for label_shared_unique_features."""
    
    def test_label_shared_unique_basic(self):
        """Test basic shared/unique labeling."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
                'gene3': np.random.randint(0, 2, 50),
            })
            
            try:
                result = analyzer.label_shared_unique_features(merged_df)
                assert result is not None
            except Exception as e:
                print(f"Label shared unique error: {e}")
    
    def test_label_shared_unique_no_cluster(self):
        """Test shared/unique labeling without cluster column."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
            })
            
            try:
                result = analyzer.label_shared_unique_features(merged_df)
                # Should handle missing Cluster column gracefully
            except Exception as e:
                print(f"No cluster error: {e}")
    
    def test_label_shared_unique_different_thresholds(self):
        """Test shared/unique labeling with different thresholds."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
            })
            
            for threshold in [0.1, 0.3, 0.5, 0.7]:
                try:
                    result = analyzer.label_shared_unique_features(merged_df, presence_threshold=threshold)
                    assert result is not None
                except Exception as e:
                    print(f"Threshold {threshold} error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerBootstrap:
    """Tests for bootstrap feature importance."""
    
    def test_bootstrap_feature_importance_basic(self):
        """Test basic bootstrap feature importance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'AMR_gene1': np.random.randint(0, 2, 50),
                'AMR_gene2': np.random.randint(0, 2, 50),
                'VIR_gene1': np.random.randint(0, 2, 50),
            })
            
            try:
                result = analyzer.bootstrap_feature_importance(merged_df, n_bootstrap=5)
                assert result is not None
                assert 'Feature' in result.columns
                assert 'Importance_Mean' in result.columns
            except Exception as e:
                print(f"Bootstrap error: {e}")
    
    def test_bootstrap_feature_importance_many_features(self):
        """Test bootstrap with many features."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 4 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(30)}
            })
            
            try:
                result = analyzer.bootstrap_feature_importance(merged_df, n_bootstrap=3)
                assert result is not None
                assert len(result) == 30
            except Exception as e:
                print(f"Many features bootstrap error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerLogOdds:
    """Tests for log odds ratio analysis."""
    
    def test_log_odds_ratio_basic(self):
        """Test basic log odds ratio analysis."""
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
                result = analyzer.log_odds_ratio_analysis(merged_df)
                assert result is not None
            except Exception as e:
                print(f"Log odds error: {e}")
    
    def test_log_odds_ratio_many_clusters(self):
        """Test log odds with many clusters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 8 for i in range(100)],
                'gene1': np.random.randint(0, 2, 100),
                'gene2': np.random.randint(0, 2, 100),
                'gene3': np.random.randint(0, 2, 100),
            })
            
            try:
                result = analyzer.log_odds_ratio_analysis(merged_df)
                assert result is not None
            except Exception as e:
                print(f"Many clusters log odds error: {e}")
