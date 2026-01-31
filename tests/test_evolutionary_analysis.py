#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for EvolutionaryAnalysis class.

Tests evolutionary metrics and phylogenetic analysis.
"""

import numpy as np
import pandas as pd
import pytest

try:
    from Bio import Phylo
    from io import StringIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

try:
    from strepsuis_phylotrait.phylo_analysis_core import EvolutionaryAnalysis
    EVO_AVAILABLE = True
except ImportError:
    EvolutionaryAnalysis = None
    EVO_AVAILABLE = False


@pytest.mark.skipif(not EVO_AVAILABLE, reason="EvolutionaryAnalysis not available")
class TestEvolutionaryAnalysis:
    """Test EvolutionaryAnalysis class."""
    
    @pytest.fixture
    def sample_tree(self):
        """Create sample phylogenetic tree."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        # Simple tree: ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);
        tree_str = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"
        tree = Phylo.read(StringIO(tree_str), "newick")
        return tree
    
    @pytest.fixture
    def evo_analysis(self):
        """Create EvolutionaryAnalysis instance."""
        return EvolutionaryAnalysis()
    
    def test_evolutionary_cluster_analysis(self, evo_analysis, sample_tree):
        """Test evolutionary cluster analysis."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        labels = np.array([1, 1, 2, 2])
        strain_names = ['A', 'B', 'C', 'D']
        mask = np.array([True, True, True, True])
        
        result = evo_analysis.evolutionary_cluster_analysis(
            sample_tree, labels, strain_names, mask
        )
        
        assert isinstance(result, (list, pd.DataFrame))
        if isinstance(result, pd.DataFrame) and not result.empty:
            assert "Cluster_ID" in result.columns
        elif isinstance(result, list) and len(result) > 0:
            assert "Cluster_ID" in result[0] or "cluster_id" in result[0]
    
    def test_calculate_phylogenetic_signal_fritz_purvis(self, evo_analysis, sample_tree):
        """Test Fritz-Purvis D statistic calculation."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        trait_data = pd.DataFrame({
            'Trait_1': [1, 1, 0, 0],
            'Trait_2': [0, 0, 1, 1],
        }, index=['A', 'B', 'C', 'D'])
        
        result = evo_analysis.calculate_phylogenetic_signal_fritz_purvis(
            sample_tree, trait_data
        )
        
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert 'trait' in result.columns
            assert 'd_statistic' in result.columns
            assert 'interpretation' in result.columns
    
    def test_sister_clade_differences(self, evo_analysis, sample_tree):
        """Test sister clade differences calculation."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        trait_vector = np.array([1, 1, 0, 0])
        strain_names = ['A', 'B', 'C', 'D']
        
        scd = evo_analysis._sister_clade_differences(
            sample_tree, strain_names, trait_vector
        )
        
        assert isinstance(scd, (int, float, np.integer, np.floating))
        assert scd >= 0
    
    def test_expected_scd_random(self, evo_analysis):
        """Test expected SCD under random distribution."""
        trait_vector = np.array([1, 0, 1, 0, 1])
        
        scd_random = evo_analysis._expected_scd_random(trait_vector)
        
        # Expected = n - 1 for random
        assert scd_random == len(trait_vector) - 1
    
    def test_expected_scd_brownian(self, evo_analysis, sample_tree):
        """Test expected SCD under Brownian motion."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        trait_vector = np.array([1, 1, 0, 0])
        strain_names = ['A', 'B', 'C', 'D']
        
        scd_brownian = evo_analysis._expected_scd_brownian(
            sample_tree, strain_names, trait_vector
        )
        
        assert isinstance(scd_brownian, (int, float, np.integer, np.floating))
        assert scd_brownian >= 0
    
    def test_calculate_beta_diversity(self, evo_analysis, sample_tree):
        """Test beta diversity calculation."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        labels = np.array([1, 1, 2, 2])
        strain_names = ['A', 'B', 'C', 'D']
        mask = np.array([True, True, True, True])
        
        beta_div = evo_analysis.calculate_beta_diversity(
            sample_tree, labels, strain_names, mask
        )
        
        assert isinstance(beta_div, pd.DataFrame)
        assert beta_div.shape[0] == beta_div.shape[1]  # Square matrix
    
    def test_calculate_evolution_rates(self, evo_analysis):
        """Test evolution rates calculation."""
        cluster_df = pd.DataFrame({
            'Cluster_ID': [1, 2, 3],
            'PD': [10.0, 15.0, 20.0],
            'InternalNodes': [5, 10, 8],
        })
        
        rates = evo_analysis.calculate_evolution_rates(cluster_df)
        
        assert isinstance(rates, pd.DataFrame)
        assert 'Cluster_ID' in rates.columns
        assert 'EvolutionRate' in rates.columns
        assert len(rates) == len(cluster_df)
    
    def test_evolutionary_cluster_analysis_single_strain(self, evo_analysis, sample_tree):
        """Test evolutionary cluster analysis with single strain per cluster."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        labels = np.array([1, 2, 3, 4])
        strain_names = ['A', 'B', 'C', 'D']
        mask = np.array([True, True, True, True])
        
        result = evo_analysis.evolutionary_cluster_analysis(
            sample_tree, labels, strain_names, mask
        )
        
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert 'Cluster_ID' in result.columns
            assert 'PD' in result.columns
