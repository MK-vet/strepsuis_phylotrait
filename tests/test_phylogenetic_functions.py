#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for phylogenetic analysis functions.

Tests phylogenetic diversity, pairwise distances, and clustering functions.
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
    from strepsuis_phylotrait.phylo_analysis_core import (
        EvolutionaryAnalysis,
    )
    EVO_AVAILABLE = True
except (ImportError, OSError):
    EvolutionaryAnalysis = None
    EVO_AVAILABLE = False


@pytest.mark.skipif(not EVO_AVAILABLE, reason="EvolutionaryAnalysis not available")
class TestPhylogeneticFunctions:
    """Test phylogenetic analysis functions."""
    
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
    def sample_trait_data(self):
        """Create sample trait data."""
        return pd.DataFrame({
            'Trait_1': [1, 1, 0, 0],
            'Trait_2': [0, 0, 1, 1],
        }, index=['A', 'B', 'C', 'D'])
    
    def test_evolutionary_cluster_analysis(self, sample_tree):
        """Test evolutionary cluster analysis."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        labels = np.array([1, 1, 2, 2])
        strain_names = ['A', 'B', 'C', 'D']
        mask = np.array([True, True, True, True])
        
        result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
            sample_tree, labels, strain_names, mask
        )
        
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert 'Cluster_ID' in result.columns
    
    def test_calculate_beta_diversity(self, sample_tree):
        """Test beta diversity calculation."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        labels = np.array([1, 1, 2, 2])
        strain_names = ['A', 'B', 'C', 'D']
        mask = np.array([True, True, True, True])
        
        result = EvolutionaryAnalysis.calculate_beta_diversity(
            sample_tree, labels, strain_names, mask
        )
        
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert result.shape[0] == result.shape[1]  # Square matrix
    
    def test_calculate_evolution_rates(self):
        """Test evolution rates calculation."""
        cluster_df = pd.DataFrame({
            'Cluster_ID': [1, 2],
            'PD': [10.0, 20.0],
            'InternalNodes': [5, 10],
        })
        
        result = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
        
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert 'Cluster_ID' in result.columns
            assert 'EvolutionRate' in result.columns
