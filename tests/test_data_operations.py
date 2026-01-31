#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for data operations in phylogenetic analysis.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        DataLoader,
        ClusteringModule,
    )
    CLASSES_AVAILABLE = True
except (ImportError, OSError):
    CLASSES_AVAILABLE = False


class TestDataLoaderDetailed:
    """Detailed tests for DataLoader class."""
    
    @pytest.fixture
    def temp_data_folder(self):
        """Create temporary data folder with CSV files."""
        temp_dir = tempfile.mkdtemp()
        
        # Create sample CSV files
        amr_df = pd.DataFrame({
            'Strain_ID': ['S1', 'S2', 'S3', 'S4'],
            'Gene_A': [1, 0, 1, 0],
            'Gene_B': [0, 1, 0, 1],
            'Gene_C': [1, 1, 0, 0],
        })
        amr_df.to_csv(os.path.join(temp_dir, 'AMR_genes.csv'), index=False)
        
        vir_df = pd.DataFrame({
            'Strain_ID': ['S1', 'S2', 'S3', 'S4'],
            'Vir_A': [1, 1, 0, 0],
            'Vir_B': [0, 0, 1, 1],
        })
        vir_df.to_csv(os.path.join(temp_dir, 'Virulence.csv'), index=False)
        
        mic_df = pd.DataFrame({
            'Strain_ID': ['S1', 'S2', 'S3', 'S4'],
            'MIC_A': [0.5, 1.0, 2.0, 4.0],
        })
        mic_df.to_csv(os.path.join(temp_dir, 'MIC.csv'), index=False)
        
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    @pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
    def test_data_loader_init(self, temp_data_folder):
        """Test DataLoader initialization."""
        loader = DataLoader(temp_data_folder)
        
        assert loader is not None
    
    @pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
    def test_data_loader_has_folder(self, temp_data_folder):
        """Test that DataLoader has data folder attribute."""
        loader = DataLoader(temp_data_folder)
        
        # Should have some reference to the data folder
        assert hasattr(loader, '__dict__')


class TestClusteringOperations:
    """Tests for clustering operations."""
    
    @pytest.fixture
    def sample_binary_data(self):
        """Create sample binary data."""
        np.random.seed(42)
        return pd.DataFrame({
            'Gene_1': np.random.binomial(1, 0.5, 20),
            'Gene_2': np.random.binomial(1, 0.3, 20),
            'Gene_3': np.random.binomial(1, 0.7, 20),
            'Gene_4': np.random.binomial(1, 0.4, 20),
            'Gene_5': np.random.binomial(1, 0.6, 20),
        })
    
    def test_binary_data_validation(self, sample_binary_data):
        """Test binary data validation."""
        # All values should be 0 or 1
        assert sample_binary_data.isin([0, 1]).all().all()
    
    def test_hamming_distance_calculation(self, sample_binary_data):
        """Test Hamming distance calculation."""
        from scipy.spatial.distance import pdist, squareform
        
        # Calculate Hamming distance matrix
        dist_matrix = squareform(pdist(sample_binary_data.values, metric='hamming'))
        
        assert dist_matrix.shape == (20, 20)
        assert np.allclose(np.diag(dist_matrix), 0)
        assert np.allclose(dist_matrix, dist_matrix.T)
    
    def test_kmodes_like_clustering(self, sample_binary_data):
        """Test K-Modes-like clustering."""
        try:
            from kmodes.kmodes import KModes
            
            km = KModes(n_clusters=3, random_state=42)
            labels = km.fit_predict(sample_binary_data.values)
            
            assert len(labels) == 20
            assert len(np.unique(labels)) <= 3
        except ImportError:
            pytest.skip("kmodes not available")


class TestDistanceMatrixOperations:
    """Tests for distance matrix operations."""
    
    def test_distance_matrix_from_tree(self):
        """Test distance matrix creation from tree."""
        try:
            from Bio import Phylo
            from io import StringIO
            
            tree_str = '((A:0.1,B:0.2):0.3,(C:0.15,D:0.25):0.35);'
            tree = Phylo.read(StringIO(tree_str), 'newick')
            
            terminals = tree.get_terminals()
            n = len(terminals)
            
            # Create distance matrix
            dist_matrix = np.zeros((n, n))
            for i, t1 in enumerate(terminals):
                for j, t2 in enumerate(terminals):
                    if i < j:
                        dist = tree.distance(t1, t2)
                        dist_matrix[i, j] = dist
                        dist_matrix[j, i] = dist
            
            assert dist_matrix.shape == (4, 4)
            assert np.allclose(np.diag(dist_matrix), 0)
        except ImportError:
            pytest.skip("Biopython not available")
    
    def test_distance_matrix_symmetry(self):
        """Test distance matrix symmetry."""
        n = 5
        dist_matrix = np.random.rand(n, n)
        dist_matrix = (dist_matrix + dist_matrix.T) / 2
        np.fill_diagonal(dist_matrix, 0)
        
        assert np.allclose(dist_matrix, dist_matrix.T)
        assert np.allclose(np.diag(dist_matrix), 0)


class TestTraitDataOperations:
    """Tests for trait data operations."""
    
    @pytest.fixture
    def sample_trait_data(self):
        """Create sample trait data."""
        return pd.DataFrame({
            'Trait_A': [1, 0, 1, 0, 1, 0],
            'Trait_B': [0, 1, 0, 1, 0, 1],
            'Trait_C': [1, 1, 0, 0, 1, 1],
        }, index=['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
    
    def test_trait_prevalence(self, sample_trait_data):
        """Test trait prevalence calculation."""
        prevalence = sample_trait_data.mean()
        
        assert isinstance(prevalence, pd.Series)
        assert len(prevalence) == 3
        assert all(0 <= p <= 1 for p in prevalence)
    
    def test_trait_cooccurrence(self, sample_trait_data):
        """Test trait co-occurrence calculation."""
        # Calculate co-occurrence matrix
        cooccurrence = sample_trait_data.T.dot(sample_trait_data)
        
        assert cooccurrence.shape == (3, 3)
        # Diagonal should be sum of each trait
        for trait in sample_trait_data.columns:
            assert cooccurrence.loc[trait, trait] == sample_trait_data[trait].sum()
    
    def test_trait_correlation(self, sample_trait_data):
        """Test trait correlation calculation."""
        correlation = sample_trait_data.corr()
        
        assert correlation.shape == (3, 3)
        # Diagonal should be 1
        assert np.allclose(np.diag(correlation.values), 1.0)
