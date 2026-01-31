#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for TreeAwareClusteringModule with real data.
"""

import os
import pytest
import pandas as pd
import numpy as np

try:
    from Bio import Phylo
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

REAL_DATA_PATH = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data')

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        TreeAwareClusteringModule,
        PhylogeneticCore,
    )
    CORE_AVAILABLE = True
except (ImportError, OSError):
    CORE_AVAILABLE = False


@pytest.fixture
def real_tree():
    """Load real SNP tree."""
    tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
    if os.path.exists(tree_path) and BIO_AVAILABLE:
        return Phylo.read(tree_path, 'newick')
    pytest.skip("Tree not found")


@pytest.fixture
def amr_data():
    """Load AMR data."""
    path = os.path.join(REAL_DATA_PATH, 'AMR_genes.csv')
    if os.path.exists(path):
        return pd.read_csv(path).set_index('Strain_ID')
    pytest.skip("AMR data not found")


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestTreeAwareClusteringReal:
    """Test TreeAwareClusteringModule with real tree."""
    
    def test_init_with_real_tree(self, real_tree):
        """Test initialization with real tree."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        assert clustering is not None
        assert clustering.tree is not None
        assert len(clustering.terminals) > 50
    
    def test_tree_cluster_algorithm_max(self, real_tree):
        """Test tree clustering with max method."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        # Create distance matrix using PhylogeneticCore
        dist_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(real_tree)
        labels = clustering.tree_cluster_algorithm(dist_matrix, method='max')
        
        assert len(labels) == len(terminals)
        assert len(np.unique(labels)) >= 1
    
    def test_tree_cluster_algorithm_sum(self, real_tree):
        """Test tree clustering with sum method."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        dist_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(real_tree)
        labels = clustering.tree_cluster_algorithm(dist_matrix, method='sum')
        
        assert len(labels) == len(terminals)
    
    def test_tree_cluster_algorithm_avg(self, real_tree):
        """Test tree clustering with avg method."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        dist_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(real_tree)
        labels = clustering.tree_cluster_algorithm(dist_matrix, method='avg')
        
        assert len(labels) == len(terminals)
    
    def test_tree_cluster_with_threshold(self, real_tree):
        """Test tree clustering with custom threshold."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        dist_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(real_tree)
        labels = clustering.tree_cluster_algorithm(dist_matrix, method='max', threshold=0.1)
        
        assert len(labels) == len(terminals)
    
    def test_auto_threshold_max(self, real_tree):
        """Test auto threshold calculation for max method."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        threshold = clustering._auto_threshold_max()
        
        assert isinstance(threshold, float)
        assert threshold > 0
    
    def test_auto_threshold_sum(self, real_tree):
        """Test auto threshold calculation for sum method."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        threshold = clustering._auto_threshold_sum()
        
        assert isinstance(threshold, float)
        assert threshold > 0
    
    def test_auto_threshold_avg(self, real_tree):
        """Test auto threshold calculation for avg method."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        threshold = clustering._auto_threshold_avg()
        
        assert isinstance(threshold, float)
        assert threshold > 0


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestTreeClusteringMethods:
    """Test clustering constraint methods."""
    
    def test_check_max_constraint(self, real_tree):
        """Test max constraint checking."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        dist_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(real_tree)
        
        # Test with first 3 terminals
        result = clustering._check_max_constraint([0, 1, 2], dist_matrix, threshold=1.0)
        
        assert result in (True, False) or isinstance(result, (bool, np.bool_))
    
    def test_check_sum_constraint(self, real_tree):
        """Test sum constraint checking."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        dist_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(real_tree)
        
        result = clustering._check_sum_constraint([0, 1, 2], dist_matrix, threshold=10.0)
        
        assert result in (True, False) or isinstance(result, (bool, np.bool_))
    
    def test_check_avg_constraint(self, real_tree):
        """Test avg constraint checking."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        dist_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(real_tree)
        
        result = clustering._check_avg_constraint([0, 1, 2], dist_matrix, threshold=1.0)
        
        assert result in (True, False) or isinstance(result, (bool, np.bool_))
    
    def test_merge_clusters(self, real_tree):
        """Test cluster merging."""
        terminals = [str(t.name) for t in real_tree.get_terminals()]
        clustering = TreeAwareClusteringModule(real_tree, terminals)
        
        n = 5
        clusters = [[i] for i in range(n)]
        cluster_map = list(range(n))
        
        clustering._merge_clusters([0, 1, 2], clusters, cluster_map)
        
        # After merge, indices 0, 1, 2 should have same cluster
        assert cluster_map[0] == cluster_map[1] == cluster_map[2]
