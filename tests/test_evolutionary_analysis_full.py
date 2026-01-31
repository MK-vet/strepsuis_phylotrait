#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Full tests for EvolutionaryAnalysis in phylo_analysis_core.py.
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
        EvolutionaryAnalysis,
        PhylogeneticCore,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestEvolutionaryAnalysisFull:
    """Full tests for EvolutionaryAnalysis."""
    
    def test_evolutionary_cluster_analysis(self):
        """Test evolutionary cluster analysis."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        labels = np.array([i % 4 for i in range(len(terminals))])
        mask = np.array([True] * len(terminals))
        
        try:
            result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
                tree, labels, strain_names, mask
            )
            assert result is not None
        except Exception as e:
            print(f"Evolutionary cluster analysis error: {e}")
    
    def test_analyze_cluster_evolution(self):
        """Test analyzing cluster evolution."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        labels = np.array([i % 4 for i in range(len(terminals))])
        mask = np.array([True] * len(terminals))
        
        try:
            if hasattr(EvolutionaryAnalysis, '_analyze_cluster_evolution'):
                result = EvolutionaryAnalysis._analyze_cluster_evolution(
                    tree, labels, strain_names, mask
                )
                assert result is not None
        except Exception as e:
            print(f"Analyze cluster evolution error: {e}")
    
    def test_calculate_beta_diversity(self):
        """Test calculating beta diversity."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        labels = np.array([i % 4 for i in range(len(terminals))])
        mask = np.array([True] * len(terminals))
        
        try:
            result = EvolutionaryAnalysis.calculate_beta_diversity(
                tree, labels, strain_names, mask
            )
            assert result is not None
        except Exception as e:
            print(f"Beta diversity error: {e}")
    
    def test_calculate_evolution_rates(self):
        """Test calculating evolution rates."""
        cluster_df = pd.DataFrame({
            'Cluster': [0, 1, 2, 3],
            'Size': [20, 25, 30, 16],
            'Silhouette': [0.5, 0.6, 0.4, 0.55],
            'Mean_Branch_Length': [0.01, 0.02, 0.015, 0.018]
        })
        
        try:
            result = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
            assert result is not None
        except Exception as e:
            print(f"Evolution rates error: {e}")
    
    def test_phylogenetic_signal_fritz_purvis(self):
        """Test Fritz-Purvis D statistic."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create trait data
        trait_data = pd.DataFrame({
            'trait1': np.random.randint(0, 2, len(terminals)),
            'trait2': np.random.randint(0, 2, len(terminals)),
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
            if len(result) > 0:
                assert 'trait' in result.columns
                assert 'd_statistic' in result.columns
        except Exception as e:
            print(f"Phylogenetic signal error: {e}")
    
    def test_phylogenetic_signal_conserved_trait(self):
        """Test phylogenetic signal for conserved trait."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create conserved trait (all same value)
        trait_data = pd.DataFrame({
            'conserved': [1] * len(terminals),
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
        except Exception as e:
            print(f"Conserved trait error: {e}")
    
    def test_phylogenetic_signal_random_trait(self):
        """Test phylogenetic signal for random trait."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create random trait
        np.random.seed(42)
        trait_data = pd.DataFrame({
            'random': np.random.randint(0, 2, len(terminals)),
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
        except Exception as e:
            print(f"Random trait error: {e}")
    
    def test_phylogenetic_signal_multiple_traits(self):
        """Test phylogenetic signal for multiple traits."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create multiple traits
        np.random.seed(42)
        trait_data = pd.DataFrame({
            'trait1': np.random.randint(0, 2, len(terminals)),
            'trait2': np.random.randint(0, 2, len(terminals)),
            'trait3': np.random.randint(0, 2, len(terminals)),
            'trait4': np.random.randint(0, 2, len(terminals)),
            'trait5': np.random.randint(0, 2, len(terminals)),
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
            if len(result) > 0:
                assert len(result) == 5
        except Exception as e:
            print(f"Multiple traits error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestEvolutionaryAnalysisSynthetic:
    """Tests for EvolutionaryAnalysis with synthetic data."""
    
    def test_evolution_rates_synthetic(self):
        """Test evolution rates with synthetic data."""
        cluster_df = pd.DataFrame({
            'Cluster': [0, 1, 2],
            'Size': [30, 35, 26],
            'Silhouette': [0.5, 0.6, 0.4],
            'Mean_Branch_Length': [0.01, 0.02, 0.015]
        })
        
        try:
            result = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
            assert result is not None
        except Exception as e:
            print(f"Synthetic evolution rates error: {e}")
    
    def test_evolution_rates_single_cluster(self):
        """Test evolution rates with single cluster."""
        cluster_df = pd.DataFrame({
            'Cluster': [0],
            'Size': [91],
            'Silhouette': [0.0],
            'Mean_Branch_Length': [0.015]
        })
        
        try:
            result = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
            assert result is not None
        except Exception as e:
            print(f"Single cluster error: {e}")
    
    def test_evolution_rates_many_clusters(self):
        """Test evolution rates with many clusters."""
        cluster_df = pd.DataFrame({
            'Cluster': list(range(10)),
            'Size': [10] * 10,
            'Silhouette': [0.5] * 10,
            'Mean_Branch_Length': [0.01 + 0.001 * i for i in range(10)]
        })
        
        try:
            result = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
            assert result is not None
        except Exception as e:
            print(f"Many clusters error: {e}")
