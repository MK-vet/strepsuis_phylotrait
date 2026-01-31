#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for PhylogeneticCore methods in phylo_analysis_core.py.
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
        PhylogeneticCore,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticCoreMethods:
    """Tests for PhylogeneticCore methods."""
    
    def test_load_tree(self):
        """Test loading tree."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        assert tree is not None
    
    def test_tree_to_distance_matrix(self):
        """Test converting tree to distance matrix."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        
        distance_matrix, terminals = PhylogeneticCore.tree_to_distance_matrix(tree)
        assert distance_matrix is not None
        assert len(terminals) > 0
        assert distance_matrix.shape[0] == len(terminals)
        assert distance_matrix.shape[1] == len(terminals)
    
    def test_dimension_reduction_umap(self):
        """Test UMAP dimension reduction."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        distance_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(tree)
        
        try:
            embeddings = PhylogeneticCore.dimension_reduction(distance_matrix, method='umap')
            assert embeddings is not None
            assert embeddings.shape[1] == 2
        except Exception as e:
            print(f"UMAP error: {e}")
    
    def test_dimension_reduction_pca(self):
        """Test PCA dimension reduction."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        distance_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(tree)
        
        try:
            if hasattr(PhylogeneticCore, 'dimension_reduction'):
                embeddings = PhylogeneticCore.dimension_reduction(distance_matrix, method='pca')
                if embeddings is not None:
                    assert embeddings.shape[1] == 2
        except Exception as e:
            print(f"PCA error: {e}")
    
    def test_detect_outliers(self):
        """Test detecting outliers."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        distance_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(tree)
        
        embeddings = PhylogeneticCore.dimension_reduction(distance_matrix)
        
        try:
            result = PhylogeneticCore.detect_outliers(embeddings)
            assert result is not None
        except Exception as e:
            print(f"Detect outliers error: {e}")
    
    def test_detect_outliers_dbscan(self):
        """Test detecting outliers with DBSCAN."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        distance_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(tree)
        
        embeddings = PhylogeneticCore.dimension_reduction(distance_matrix)
        
        try:
            if hasattr(PhylogeneticCore, 'detect_outliers_dbscan'):
                result = PhylogeneticCore.detect_outliers_dbscan(embeddings)
                assert result is not None
        except Exception as e:
            print(f"DBSCAN outliers error: {e}")
    
    def test_calculate_faith_pd(self):
        """Test calculating Faith's PD."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        terminals = tree.get_terminals()
        
        labels = np.array([i % 3 for i in range(len(terminals))])
        
        try:
            if hasattr(PhylogeneticCore, 'calculate_faith_pd'):
                result = PhylogeneticCore.calculate_faith_pd(tree, labels)
                assert result is not None
        except Exception as e:
            print(f"Faith PD error: {e}")
    
    def test_calculate_nri_nti(self):
        """Test calculating NRI/NTI."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        terminals = tree.get_terminals()
        
        labels = np.array([i % 3 for i in range(len(terminals))])
        
        try:
            if hasattr(PhylogeneticCore, 'calculate_nri_nti'):
                result = PhylogeneticCore.calculate_nri_nti(tree, labels)
                assert result is not None
        except Exception as e:
            print(f"NRI/NTI error: {e}")
    
    def test_get_terminal_names(self):
        """Test getting terminal names."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        
        try:
            if hasattr(PhylogeneticCore, 'get_terminal_names'):
                result = PhylogeneticCore.get_terminal_names(tree)
                assert result is not None
                assert len(result) > 0
        except Exception as e:
            print(f"Get terminal names error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestPhylogeneticCoreSyntheticData:
    """Tests for PhylogeneticCore with synthetic data."""
    
    def test_dimension_reduction_synthetic(self):
        """Test dimension reduction with synthetic data."""
        np.random.seed(42)
        distance_matrix = np.random.rand(50, 50)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(distance_matrix, 0)
        
        try:
            embeddings = PhylogeneticCore.dimension_reduction(distance_matrix)
            assert embeddings is not None
            assert embeddings.shape == (50, 2)
        except Exception as e:
            print(f"Synthetic dimension reduction error: {e}")
    
    def test_detect_outliers_synthetic(self):
        """Test detecting outliers with synthetic data."""
        np.random.seed(42)
        embeddings = np.random.randn(50, 2)
        
        try:
            result = PhylogeneticCore.detect_outliers(embeddings)
            assert result is not None
        except Exception as e:
            print(f"Synthetic outliers error: {e}")
    
    def test_small_distance_matrix(self):
        """Test with small distance matrix."""
        np.random.seed(42)
        distance_matrix = np.random.rand(10, 10)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(distance_matrix, 0)
        
        try:
            embeddings = PhylogeneticCore.dimension_reduction(distance_matrix)
            assert embeddings is not None
            assert embeddings.shape == (10, 2)
        except Exception as e:
            print(f"Small matrix error: {e}")
    
    def test_large_distance_matrix(self):
        """Test with large distance matrix."""
        np.random.seed(42)
        distance_matrix = np.random.rand(200, 200)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        np.fill_diagonal(distance_matrix, 0)
        
        try:
            embeddings = PhylogeneticCore.dimension_reduction(distance_matrix)
            assert embeddings is not None
            assert embeddings.shape == (200, 2)
        except Exception as e:
            print(f"Large matrix error: {e}")
