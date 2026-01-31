"""
Tests for phylo_analysis_core module.

Tests main classes and functions for phylogenetic analysis.
"""

import numpy as np
import pandas as pd
import pytest
import tempfile
import os
import logging

from strepsuis_phylotrait.phylo_analysis_core import (
    print_memory_usage,
    print_section_header,
    print_step,
    create_template_directory,
    ParallelProcessor,
    PhylogeneticCore,
    ClusteringModule,
    EvolutionaryAnalysis,
    DataLoader,
    Visualizer,
    TraitAnalyzer,
    MCAAnalyzer,
)


class TestUtilityFunctions:
    """Tests for utility functions."""
    
    def test_print_memory_usage(self):
        """Test memory usage printing."""
        mem = print_memory_usage()
        
        assert mem > 0
    
    def test_print_section_header(self):
        """Test section header printing."""
        # Should not raise
        print_section_header("Test Section")
    
    def test_print_step(self):
        """Test step printing."""
        # Should not raise
        print_step(1, 5, "Test step")


class TestCreateTemplateDirectory:
    """Tests for template directory creation."""
    
    def test_create_template_directory(self):
        """Test template directory creation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            old_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                create_template_directory()
                
                assert os.path.exists("templates")
            finally:
                os.chdir(old_cwd)


class TestParallelProcessor:
    """Tests for ParallelProcessor class."""
    
    def test_parallel_feature_importance(self):
        """Test parallel feature importance."""
        np.random.seed(42)
        X = np.random.randint(0, 2, (50, 5))
        y = np.random.randint(0, 2, 50)
        
        importances = ParallelProcessor.parallel_feature_importance(
            X, y, n_bootstrap=3, n_jobs=1
        )
        
        assert importances.shape == (3, 5)


class TestPhylogeneticCore:
    """Tests for PhylogeneticCore class."""
    
    def test_dimension_reduction(self):
        """Test UMAP dimension reduction."""
        np.random.seed(42)
        # Create a symmetric distance matrix
        n = 20
        matrix = np.random.rand(n, n)
        matrix = (matrix + matrix.T) / 2
        np.fill_diagonal(matrix, 0)
        
        embeddings = PhylogeneticCore.dimension_reduction(
            matrix, n_components=2, n_neighbors=5
        )
        
        assert embeddings.shape == (n, 2)
    
    def test_detect_outliers(self):
        """Test outlier detection."""
        np.random.seed(42)
        embeddings = np.random.rand(50, 2)
        
        filtered, mask = PhylogeneticCore.detect_outliers(
            embeddings, contamination=0.1
        )
        
        assert len(mask) == 50
        assert len(filtered) <= 50


class TestClusteringModule:
    """Tests for ClusteringModule class."""
    
    def test_initialization(self):
        """Test ClusteringModule initialization."""
        module = ClusteringModule(n_clusters_range=(2, 5), seed=42)
        
        assert module.n_clusters_range == (2, 5)
        assert module.seed == 42
    
    def test_assign_outliers_to_clusters(self):
        """Test outlier assignment."""
        np.random.seed(42)
        embeddings = np.random.rand(20, 2)
        mask = np.array([True] * 18 + [False, False])
        labels = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        
        module = ClusteringModule()
        assignments = module.assign_outliers_to_clusters(embeddings, mask, labels)
        
        assert len(assignments) == 2  # 2 outliers


class TestEvolutionaryAnalysis:
    """Tests for EvolutionaryAnalysis class."""
    
    def test_static_methods_exist(self):
        """Test that static methods exist."""
        assert hasattr(EvolutionaryAnalysis, 'evolutionary_cluster_analysis')


class TestDataLoader:
    """Tests for DataLoader class."""
    
    def test_initialization(self):
        """Test DataLoader initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = DataLoader(tmpdir)
            
            assert loader.data_folder == tmpdir


class TestTraitAnalyzer:
    """Tests for TraitAnalyzer class."""
    
    def test_initialization(self):
        """Test TraitAnalyzer initialization."""
        np.random.seed(42)
        data = pd.DataFrame({
            'Strain_ID': [f'S{i}' for i in range(20)],
            'Trait1': np.random.randint(0, 2, 20),
            'Trait2': np.random.randint(0, 2, 20)
        })
        
        analyzer = TraitAnalyzer(data, output_folder='.')
        
        assert analyzer is not None


class TestMCAAnalyzer:
    """Tests for MCAAnalyzer class."""
    
    def test_initialization(self):
        """Test MCAAnalyzer initialization."""
        np.random.seed(42)
        data = pd.DataFrame(np.random.randint(0, 2, (30, 5)))
        
        analyzer = MCAAnalyzer(data, output_folder='.')
        
        assert analyzer is not None


class TestVisualizer:
    """Tests for Visualizer class."""
    
    def test_initialization(self):
        """Test Visualizer initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(output_folder=tmpdir)
            
            assert viz.output_folder == tmpdir


class TestIntegration:
    """Integration tests."""
    
    def test_clustering_module_creation(self):
        """Test clustering module creation."""
        module = ClusteringModule(n_clusters_range=(2, 4), seed=42)
        
        assert module.n_clusters_range == (2, 4)
    
    def test_phylogenetic_core_methods(self):
        """Test PhylogeneticCore methods."""
        np.random.seed(42)
        n = 15
        matrix = np.random.rand(n, n)
        matrix = (matrix + matrix.T) / 2
        np.fill_diagonal(matrix, 0)
        
        # Test dimension reduction
        embeddings = PhylogeneticCore.dimension_reduction(
            matrix, n_components=2, n_neighbors=5
        )
        
        assert embeddings.shape == (n, 2)
        
        # Test outlier detection
        filtered, mask = PhylogeneticCore.detect_outliers(embeddings)
        
        assert len(mask) == n
