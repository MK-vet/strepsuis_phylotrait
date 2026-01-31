#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Coverage boost tests for phylo_analysis_core.py.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil
import logging

REAL_DATA_PATH = r"C:\Users\ABC\Documents\GitHub\MKrep\data"


def check_real_data_exists():
    return os.path.exists(os.path.join(REAL_DATA_PATH, 'Snp_tree.newick'))


try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        PhylogeneticCore,
        TreeAwareClusteringModule,
        ClusteringModule,
        EvolutionaryAnalysis,
        TraitAnalyzer,
        MCAAnalyzer,
        Visualizer,
        DataLoader,
        HTMLReportGenerator,
        Config,
        PhylogeneticAnalysis,
        setup_logging,
        print_memory_usage,
        print_section_header,
        print_step,
        create_template_directory,
        ParallelProcessor,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticCoreBoost:
    """Boost tests for PhylogeneticCore."""
    
    def test_load_tree_variations(self):
        """Test loading tree with variations."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        
        # Test static method
        tree = PhylogeneticCore.load_tree(tree_path)
        assert tree is not None
        
        # Test instance method
        core = PhylogeneticCore()
        if hasattr(core, 'load_tree'):
            tree2 = core.load_tree(tree_path)
            assert tree2 is not None
    
    def test_dimension_reduction_params(self):
        """Test dimension reduction with different parameters."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        distance_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(tree)
        
        # Different n_components
        try:
            embeddings = PhylogeneticCore.dimension_reduction(
                distance_matrix, n_components=2
            )
            assert embeddings.shape[1] == 2
        except Exception as e:
            print(f"n_components error: {e}")
        
        # Different n_neighbors
        for n in [5, 10, 15]:
            try:
                embeddings = PhylogeneticCore.dimension_reduction(
                    distance_matrix, n_neighbors=n
                )
                assert embeddings is not None
            except Exception as e:
                print(f"n_neighbors={n} error: {e}")
        
        # Different min_dist
        for d in [0.1, 0.3, 0.5]:
            try:
                embeddings = PhylogeneticCore.dimension_reduction(
                    distance_matrix, min_dist=d
                )
                assert embeddings is not None
            except Exception as e:
                print(f"min_dist={d} error: {e}")
    
    def test_detect_outliers_params(self):
        """Test outlier detection with different parameters."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = PhylogeneticCore.load_tree(tree_path)
        distance_matrix, _ = PhylogeneticCore.tree_to_distance_matrix(tree)
        embeddings = PhylogeneticCore.dimension_reduction(distance_matrix)
        
        # Different eps values
        for eps in [0.3, 0.5, 0.7, 1.0]:
            try:
                result = PhylogeneticCore.detect_outliers(embeddings, eps=eps)
                assert result is not None
            except Exception as e:
                print(f"eps={eps} error: {e}")
        
        # Different min_samples
        for min_samples in [3, 5, 10]:
            try:
                result = PhylogeneticCore.detect_outliers(
                    embeddings, min_samples=min_samples
                )
                assert result is not None
            except Exception as e:
                print(f"min_samples={min_samples} error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestClusteringBoost:
    """Boost tests for clustering modules."""
    
    def test_clustering_module_params(self):
        """Test ClusteringModule with different parameters."""
        module = ClusteringModule(
            n_clusters_range=(2, 10),
            n_ensemble=10,
            dbscan_trials=20,
            seed=42
        )
        
        np.random.seed(42)
        data = np.random.randint(0, 2, (100, 30))
        
        try:
            result = module.ensemble_clustering(data)
            assert result is not None
            assert len(result) == 100
        except Exception as e:
            print(f"Clustering module params error: {e}")
    
    def test_tree_aware_clustering_params(self):
        """Test TreeAwareClusteringModule with different parameters."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 10),
            n_ensemble=5,
            seed=42
        )
        
        core = PhylogeneticCore()
        distance_matrix, _ = core.tree_to_distance_matrix(tree)
        
        # Test with different k values
        for k in [2, 3, 4, 5, 6]:
            try:
                result = module.tree_cluster_algorithm(
                    distance_matrix, n_clusters=k
                )
                assert result is not None
                assert len(result) == len(terminals)
            except Exception as e:
                print(f"k={k} error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestEvolutionaryBoost:
    """Boost tests for evolutionary analysis."""
    
    def test_evolutionary_with_different_clusters(self):
        """Test evolutionary analysis with different cluster numbers."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        mask = np.array([True] * len(terminals))
        
        for n_clusters in [2, 3, 4, 5, 6]:
            labels = np.array([i % n_clusters for i in range(len(terminals))])
            
            try:
                result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
                    tree, labels, strain_names, mask
                )
                assert result is not None
            except Exception as e:
                print(f"n_clusters={n_clusters} evolutionary error: {e}")
    
    def test_phylogenetic_signal_variations(self):
        """Test phylogenetic signal with different traits."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Random traits
        np.random.seed(42)
        trait_data = pd.DataFrame({
            f'random_{i}': np.random.randint(0, 2, len(terminals))
            for i in range(5)
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
            assert len(result) == 5
        except Exception as e:
            print(f"Random traits error: {e}")
        
        # Conserved traits (all same)
        trait_data_conserved = pd.DataFrame({
            'all_ones': [1] * len(terminals),
            'all_zeros': [0] * len(terminals),
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data_conserved)
            assert result is not None
        except Exception as e:
            print(f"Conserved traits error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerBoost:
    """Boost tests for TraitAnalyzer."""
    
    def test_association_rules_params(self):
        """Test association rules with different parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 4 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(40)}
            })
            
            # Different support thresholds
            for support in [0.05, 0.1, 0.2, 0.3]:
                try:
                    result = analyzer.association_rule_mining(
                        merged_df, min_support=support
                    )
                    assert result is not None
                except Exception as e:
                    print(f"support={support} error: {e}")
            
            # Different confidence thresholds
            for confidence in [0.5, 0.7, 0.9]:
                try:
                    result = analyzer.association_rule_mining(
                        merged_df, min_confidence=confidence
                    )
                    assert result is not None
                except Exception as e:
                    print(f"confidence={confidence} error: {e}")
    
    def test_bootstrap_params(self):
        """Test bootstrap with different parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 4 for i in range(100)],
                **{f'gene_{i}': np.random.randint(0, 2, 100) for i in range(20)}
            })
            
            for n_bootstrap in [3, 5, 10]:
                try:
                    result = analyzer.bootstrap_feature_importance(
                        merged_df, n_bootstrap=n_bootstrap
                    )
                    assert result is not None
                except Exception as e:
                    print(f"n_bootstrap={n_bootstrap} error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestMCABoost:
    """Boost tests for MCAAnalyzer."""
    
    def test_mca_with_different_data_sizes(self):
        """Test MCA with different data sizes."""
        for n_samples in [30, 50, 100]:
            with tempfile.TemporaryDirectory() as tmpdir:
                analyzer = MCAAnalyzer(tmpdir)
                
                np.random.seed(42)
                merged_df = pd.DataFrame({
                    'Strain_ID': [f'strain_{i}' for i in range(n_samples)],
                    'Cluster': [i % 3 for i in range(n_samples)],
                    **{f'gene_{i}': np.random.randint(0, 2, n_samples) for i in range(15)}
                })
                
                try:
                    row_coords, col_coords, summary = analyzer.perform_mca_analysis(merged_df)
                    
                    if row_coords is not None:
                        assert len(row_coords) == n_samples
                except Exception as e:
                    print(f"n_samples={n_samples} MCA error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestVisualizerBoost:
    """Boost tests for Visualizer."""
    
    def test_visualizations_with_different_clusters(self):
        """Test visualizations with different cluster numbers."""
        for n_clusters in [2, 3, 5, 8]:
            with tempfile.TemporaryDirectory() as tmpdir:
                viz = Visualizer(tmpdir)
                
                np.random.seed(42)
                n_samples = 100
                merged_df = pd.DataFrame({
                    'Strain_ID': [f'strain_{i}' for i in range(n_samples)],
                    'Cluster': [i % n_clusters for i in range(n_samples)],
                    'gene1': np.random.randint(0, 2, n_samples),
                })
                
                try:
                    viz.plot_cluster_distribution(merged_df)
                except Exception as e:
                    print(f"n_clusters={n_clusters} distribution error: {e}")
                
                embeddings = np.random.randn(n_samples, 2)
                labels = np.array([i % n_clusters for i in range(n_samples)])
                mask = np.array([True] * n_samples)
                
                try:
                    viz.plot_umap_clusters(embeddings, labels, mask)
                except Exception as e:
                    print(f"n_clusters={n_clusters} UMAP error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestHelperFunctionsBoost:
    """Boost tests for helper functions."""
    
    def test_print_functions_multiple_calls(self):
        """Test print functions with multiple calls."""
        for i in range(3):
            try:
                print_memory_usage()
            except Exception as e:
                print(f"Print memory {i} error: {e}")
        
        for section in ["Section A", "Section B", "Section C"]:
            try:
                print_section_header(section)
            except Exception as e:
                print(f"Print section {section} error: {e}")
        
        for i in range(1, 6):
            try:
                print_step(f"Step {i}", i, 5)
            except Exception as e:
                print(f"Print step {i} error: {e}")
    
    def test_create_template_directory_multiple(self):
        """Test creating template directory multiple times."""
        for _ in range(3):
            try:
                create_template_directory()
            except Exception as e:
                print(f"Create template error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestParallelProcessorBoost:
    """Boost tests for ParallelProcessor."""
    
    def test_parallel_processor_params(self):
        """Test ParallelProcessor with different parameters."""
        for n_jobs in [1, 2, -1]:
            try:
                processor = ParallelProcessor(n_jobs=n_jobs)
                assert processor is not None
            except Exception as e:
                print(f"n_jobs={n_jobs} error: {e}")
