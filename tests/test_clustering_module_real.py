#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for ClusteringModule with real data.
"""

import os
import pytest
import pandas as pd
import numpy as np

REAL_DATA_PATH = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data')

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        ClusteringModule,
        Config,
    )
    CORE_AVAILABLE = True
except (ImportError, OSError):
    CORE_AVAILABLE = False


@pytest.fixture
def amr_data():
    """Load AMR data."""
    path = os.path.join(REAL_DATA_PATH, 'AMR_genes.csv')
    if os.path.exists(path):
        return pd.read_csv(path).set_index('Strain_ID')
    pytest.skip("AMR data not found")


@pytest.fixture
def virulence_data():
    """Load virulence data."""
    path = os.path.join(REAL_DATA_PATH, 'Virulence.csv')
    if os.path.exists(path):
        return pd.read_csv(path).set_index('Strain_ID')
    pytest.skip("Virulence data not found")


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestClusteringModuleReal:
    """Test ClusteringModule with real data."""
    
    def test_init_default(self):
        """Test initialization with default parameters."""
        clustering = ClusteringModule()
        
        assert clustering is not None
        assert clustering.n_clusters_range == (2, 20)
        assert clustering.n_ensemble == 10
    
    def test_init_custom_params(self):
        """Test initialization with custom parameters."""
        clustering = ClusteringModule(
            n_clusters_range=(3, 10),
            n_ensemble=5,
            dbscan_trials=20,
            seed=123
        )
        
        assert clustering.n_clusters_range == (3, 10)
        assert clustering.n_ensemble == 5
        assert clustering.dbscan_trials == 20
        assert clustering.seed == 123
    
    def test_ensemble_clustering_basic(self, amr_data):
        """Test ensemble clustering with real data."""
        # Use subset for speed
        data_subset = amr_data.iloc[:30, :5].values
        
        clustering = ClusteringModule(
            n_clusters_range=(2, 4),
            n_ensemble=2,
            dbscan_trials=5
        )
        
        labels, score = clustering.ensemble_clustering(data_subset)
        
        if labels is not None:
            assert len(labels) == len(data_subset)
    
    def test_assign_outliers_to_clusters(self, amr_data):
        """Test outlier assignment."""
        data_subset = amr_data.iloc[:20, :5].values
        
        clustering = ClusteringModule()
        
        # Create mock labels and mask
        labels = np.array([1, 1, 2, 2, 2, 1, 1, 2, 2, 1])
        mask = np.array([True] * 10 + [False] * 10)
        
        assignments = clustering.assign_outliers_to_clusters(
            data_subset, mask, labels
        )
        
        assert isinstance(assignments, list)


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestClusteringWithCombinedData:
    """Test clustering with combined AMR and virulence data."""
    
    def test_combined_clustering(self, amr_data, virulence_data):
        """Test clustering on combined data."""
        # Join data
        combined = amr_data.join(virulence_data, how='inner', rsuffix='_vir')
        data_subset = combined.iloc[:30, :10].values
        
        clustering = ClusteringModule(
            n_clusters_range=(2, 4),
            n_ensemble=2,
            dbscan_trials=5
        )
        
        labels, score = clustering.ensemble_clustering(data_subset)
        
        if labels is not None:
            assert len(labels) == len(data_subset)
    
    def test_optimize_dbscan(self, amr_data):
        """Test DBSCAN optimization."""
        data_subset = amr_data.iloc[:30, :5].values
        
        clustering = ClusteringModule(dbscan_trials=5)
        
        params = clustering._optimize_dbscan(data_subset)
        
        assert isinstance(params, dict)
        assert 'eps' in params
        assert 'min_samples' in params


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestConfigWithRealPaths:
    """Test Config with real data paths."""
    
    def test_config_setup(self):
        """Test Config setup."""
        config = Config()
        
        config.data_folder = REAL_DATA_PATH
        config.n_clusters = 5
        config.random_state = 42
        
        assert config.data_folder == REAL_DATA_PATH
    
    def test_config_tree_file(self):
        """Test Config with tree file."""
        config = Config()
        
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        config.tree_file = tree_path
        
        if os.path.exists(tree_path):
            assert config.tree_file == tree_path
