#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Basic tests that don't require full module import.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np


class TestBasicDataOperations:
    """Test basic data operations without full module import."""
    
    def test_dataframe_creation(self):
        """Test basic DataFrame creation."""
        df = pd.DataFrame({
            'Strain_ID': ['S1', 'S2', 'S3'],
            'Gene_A': [1, 0, 1],
            'Gene_B': [0, 1, 0],
        })
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert 'Strain_ID' in df.columns
    
    def test_binary_data_validation(self):
        """Test binary data validation."""
        df = pd.DataFrame({
            'Gene_A': [1, 0, 1, 0],
            'Gene_B': [0, 1, 0, 1],
        })
        
        # Check all values are 0 or 1
        assert df.isin([0, 1]).all().all()
    
    def test_distance_matrix_creation(self):
        """Test distance matrix creation."""
        n = 4
        dist_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(n):
                if i != j:
                    dist_matrix[i, j] = abs(i - j) * 0.1
        
        assert dist_matrix.shape == (n, n)
        assert np.allclose(np.diag(dist_matrix), 0)
        assert np.allclose(dist_matrix, dist_matrix.T)
    
    def test_cluster_labels_creation(self):
        """Test cluster labels creation."""
        n_samples = 10
        n_clusters = 3
        
        labels = np.random.randint(0, n_clusters, n_samples)
        
        assert len(labels) == n_samples
        assert len(np.unique(labels)) <= n_clusters


class TestFileOperations:
    """Test file operations."""
    
    @pytest.fixture
    def temp_folder(self):
        """Create temporary folder."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_csv_save_and_load(self, temp_folder):
        """Test CSV save and load."""
        df = pd.DataFrame({
            'A': [1, 2, 3],
            'B': [4, 5, 6],
        })
        
        csv_path = os.path.join(temp_folder, 'test.csv')
        df.to_csv(csv_path, index=False)
        
        loaded_df = pd.read_csv(csv_path)
        
        assert len(loaded_df) == len(df)
        assert list(loaded_df.columns) == list(df.columns)
    
    def test_folder_creation(self, temp_folder):
        """Test folder creation."""
        new_folder = os.path.join(temp_folder, 'new_subfolder')
        os.makedirs(new_folder, exist_ok=True)
        
        assert os.path.exists(new_folder)
    
    def test_newick_file_creation(self, temp_folder):
        """Test Newick tree file creation."""
        tree_path = os.path.join(temp_folder, 'test.newick')
        tree_content = '((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);'
        
        with open(tree_path, 'w') as f:
            f.write(tree_content)
        
        assert os.path.exists(tree_path)
        
        with open(tree_path, 'r') as f:
            content = f.read()
        
        assert content == tree_content


class TestStatisticalOperations:
    """Test statistical operations."""
    
    def test_mean_calculation(self):
        """Test mean calculation."""
        data = np.array([1, 2, 3, 4, 5])
        
        mean = np.mean(data)
        
        assert mean == 3.0
    
    def test_std_calculation(self):
        """Test standard deviation calculation."""
        data = np.array([1, 2, 3, 4, 5])
        
        std = np.std(data)
        
        assert std > 0
    
    def test_correlation_calculation(self):
        """Test correlation calculation."""
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 6, 8, 10])
        
        corr = np.corrcoef(x, y)[0, 1]
        
        assert np.isclose(corr, 1.0)
    
    def test_percentile_calculation(self):
        """Test percentile calculation."""
        data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        
        p50 = np.percentile(data, 50)
        
        assert p50 == 5.5


class TestClusteringOperations:
    """Test clustering-related operations."""
    
    def test_hamming_distance(self):
        """Test Hamming distance calculation."""
        x = np.array([1, 0, 1, 0])
        y = np.array([1, 1, 0, 0])
        
        # Hamming distance = number of different positions / total positions
        hamming = np.sum(x != y) / len(x)
        
        assert hamming == 0.5
    
    def test_silhouette_like_metric(self):
        """Test silhouette-like metric calculation."""
        # Simple example with 2 clusters
        data = np.array([
            [0, 0],
            [0.1, 0.1],
            [10, 10],
            [10.1, 10.1],
        ])
        labels = np.array([0, 0, 1, 1])
        
        # Calculate intra-cluster distances
        cluster_0 = data[labels == 0]
        cluster_1 = data[labels == 1]
        
        intra_0 = np.mean(np.linalg.norm(cluster_0[0] - cluster_0[1]))
        intra_1 = np.mean(np.linalg.norm(cluster_1[0] - cluster_1[1]))
        
        assert intra_0 < 1  # Points in cluster 0 are close
        assert intra_1 < 1  # Points in cluster 1 are close


class TestBioPhyloCompatibility:
    """Test Bio.Phylo compatibility."""
    
    def test_biopython_import(self):
        """Test if Biopython can be imported."""
        try:
            from Bio import Phylo
            from io import StringIO
            
            tree_str = '((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);'
            tree = Phylo.read(StringIO(tree_str), 'newick')
            
            terminals = tree.get_terminals()
            
            assert len(terminals) == 4
        except ImportError:
            pytest.skip("Biopython not available")
    
    def test_tree_structure(self):
        """Test tree structure parsing."""
        try:
            from Bio import Phylo
            from io import StringIO
            
            tree_str = '((A:0.1,B:0.2):0.3,(C:0.15,D:0.25):0.35);'
            tree = Phylo.read(StringIO(tree_str), 'newick')
            
            # Get terminal names
            terminal_names = [t.name for t in tree.get_terminals()]
            
            assert 'A' in terminal_names
            assert 'B' in terminal_names
            assert 'C' in terminal_names
            assert 'D' in terminal_names
        except ImportError:
            pytest.skip("Biopython not available")
