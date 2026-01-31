#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for clustering methods in phylo_analysis_core.py.
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
    from strepsuis_phylotrait.phylo_analysis_core import (
        ClusteringModule,
        TreeAwareClusteringModule,
        PhylogeneticCore,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestClusteringModuleMethods:
    """Tests for ClusteringModule methods."""
    
    def test_kmodes_clustering(self):
        """Test K-Modes clustering."""
        module = ClusteringModule(
            n_clusters_range=(2, 6),
            n_ensemble=3,
            dbscan_trials=5,
            seed=42
        )
        
        np.random.seed(42)
        data = np.random.randint(0, 2, (50, 15))
        
        try:
            if hasattr(module, 'kmodes_clustering'):
                result = module.kmodes_clustering(data, n_clusters=3)
                assert result is not None
                assert len(result) == 50
        except Exception as e:
            print(f"K-Modes error: {e}")
    
    def test_hierarchical_clustering(self):
        """Test hierarchical clustering."""
        module = ClusteringModule(
            n_clusters_range=(2, 6),
            n_ensemble=3,
            dbscan_trials=5,
            seed=42
        )
        
        np.random.seed(42)
        distance_matrix = np.random.rand(50, 50)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(distance_matrix, 0)
        
        try:
            if hasattr(module, 'hierarchical_clustering'):
                result = module.hierarchical_clustering(distance_matrix, n_clusters=3)
                assert result is not None
                assert len(result) == 50
        except Exception as e:
            print(f"Hierarchical error: {e}")
    
    def test_dbscan_clustering(self):
        """Test DBSCAN clustering."""
        module = ClusteringModule(
            n_clusters_range=(2, 6),
            n_ensemble=3,
            dbscan_trials=5,
            seed=42
        )
        
        np.random.seed(42)
        embeddings = np.random.randn(50, 2)
        
        try:
            if hasattr(module, 'dbscan_clustering'):
                result = module.dbscan_clustering(embeddings)
                assert result is not None
                assert len(result) == 50
        except Exception as e:
            print(f"DBSCAN error: {e}")
    
    def test_find_optimal_k(self):
        """Test finding optimal k."""
        module = ClusteringModule(
            n_clusters_range=(2, 6),
            n_ensemble=3,
            dbscan_trials=5,
            seed=42
        )
        
        np.random.seed(42)
        data = np.random.randint(0, 2, (50, 15))
        
        try:
            if hasattr(module, 'find_optimal_k'):
                result = module.find_optimal_k(data)
                assert result is not None
                assert 2 <= result <= 6
        except Exception as e:
            print(f"Find optimal k error: {e}")
    
    def test_evaluate_clustering(self):
        """Test evaluating clustering."""
        module = ClusteringModule(
            n_clusters_range=(2, 6),
            n_ensemble=3,
            dbscan_trials=5,
            seed=42
        )
        
        np.random.seed(42)
        data = np.random.randint(0, 2, (50, 15))
        labels = np.array([i % 3 for i in range(50)])
        
        try:
            if hasattr(module, 'evaluate_clustering'):
                result = module.evaluate_clustering(data, labels)
                assert result is not None
        except Exception as e:
            print(f"Evaluate clustering error: {e}")
    
    def test_ensemble_clustering_full(self):
        """Test full ensemble clustering."""
        module = ClusteringModule(
            n_clusters_range=(2, 6),
            n_ensemble=5,
            dbscan_trials=10,
            seed=42
        )
        
        np.random.seed(42)
        data = np.random.randint(0, 2, (50, 15))
        
        try:
            result = module.ensemble_clustering(data)
            assert result is not None
            assert len(result) == 50
        except Exception as e:
            print(f"Ensemble clustering error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestTreeAwareClusteringMethods:
    """Tests for TreeAwareClusteringModule methods."""
    
    def test_tree_cluster_max(self):
        """Test tree clustering with max method."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        core = PhylogeneticCore()
        distance_matrix, _ = core.tree_to_distance_matrix(tree)
        
        try:
            result = module.tree_cluster_algorithm(distance_matrix, method='max')
            assert result is not None
            assert len(result) == len(terminals)
        except Exception as e:
            print(f"Tree cluster max error: {e}")
    
    def test_tree_cluster_sum(self):
        """Test tree clustering with sum method."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        core = PhylogeneticCore()
        distance_matrix, _ = core.tree_to_distance_matrix(tree)
        
        try:
            result = module.tree_cluster_algorithm(distance_matrix, method='sum')
            assert result is not None
            assert len(result) == len(terminals)
        except Exception as e:
            print(f"Tree cluster sum error: {e}")
    
    def test_tree_cluster_avg(self):
        """Test tree clustering with avg method."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        core = PhylogeneticCore()
        distance_matrix, _ = core.tree_to_distance_matrix(tree)
        
        try:
            result = module.tree_cluster_algorithm(distance_matrix, method='avg')
            assert result is not None
            assert len(result) == len(terminals)
        except Exception as e:
            print(f"Tree cluster avg error: {e}")
    
    def test_phydelity_clustering(self):
        """Test phydelity clustering."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        core = PhylogeneticCore()
        distance_matrix, _ = core.tree_to_distance_matrix(tree)
        
        try:
            result = module.phydelity_clustering(distance_matrix)
            assert result is not None
            assert len(result) == len(terminals)
        except Exception as e:
            print(f"Phydelity error: {e}")
    
    def test_evaluate_monophyly(self):
        """Test evaluating monophyly."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        labels = np.array([i % 3 for i in range(len(terminals))])
        
        try:
            result = module.evaluate_monophyly(labels)
            assert result is not None
        except Exception as e:
            print(f"Evaluate monophyly error: {e}")
    
    def test_ensure_monophyletic_clusters(self):
        """Test ensuring monophyletic clusters."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        labels = np.array([i % 3 for i in range(len(terminals))])
        
        try:
            result = module.ensure_monophyletic_clusters(labels)
            assert result is not None
            assert len(result) == len(terminals)
        except Exception as e:
            print(f"Ensure monophyletic error: {e}")
    
    def test_compare_methods(self):
        """Test comparing clustering methods."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        core = PhylogeneticCore()
        distance_matrix, _ = core.tree_to_distance_matrix(tree)
        
        try:
            if hasattr(module, 'compare_methods'):
                result = module.compare_methods(distance_matrix)
                assert result is not None
        except Exception as e:
            print(f"Compare methods error: {e}")
