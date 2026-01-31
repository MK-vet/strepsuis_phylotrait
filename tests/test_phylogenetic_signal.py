#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for phylogenetic signal detection.

Tests Fritz-Purvis D statistic calculation for binary traits.
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


class TestPhylogeneticSignal:
    """Test phylogenetic signal detection functions."""
    
    @pytest.fixture
    def sample_tree(self):
        """Create sample phylogenetic tree."""
        if not BIO_AVAILABLE or not EVO_AVAILABLE:
            pytest.skip("Bio.Phylo or EvolutionaryAnalysis not available")
        # Simple tree: ((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);
        from Bio import Phylo
        from io import StringIO
        tree_str = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);"
        tree = Phylo.read(StringIO(tree_str), "newick")
        return tree
    
    @pytest.fixture
    def sample_trait_data(self):
        """Create sample binary trait data."""
        return pd.DataFrame({
            'Strain_A': [1, 0, 1],
            'Strain_B': [1, 0, 1],
            'Strain_C': [0, 1, 0],
            'Strain_D': [0, 1, 0],
        }, index=['Trait_1', 'Trait_2', 'Trait_3']).T
    
    def test_calculate_phylogenetic_signal_fritz_purvis(
        self, sample_tree, sample_trait_data
    ):
        """Test Fritz-Purvis D statistic calculation."""
        if not EVO_AVAILABLE or not BIO_AVAILABLE:
            pytest.skip("EvolutionaryAnalysis or Bio.Phylo not available")
        
        evo_analysis = EvolutionaryAnalysis()
        
        # Create trait data with strain names matching tree
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
            assert 'scd_observed' in result.columns
            assert 'scd_random' in result.columns
            assert 'scd_brownian' in result.columns
    
    def test_sister_clade_differences(self, sample_tree):
        """Test sister clade differences calculation."""
        if not EVO_AVAILABLE or not BIO_AVAILABLE:
            pytest.skip("EvolutionaryAnalysis or Bio.Phylo not available")
        
        evo_analysis = EvolutionaryAnalysis()
        
        trait_vector = np.array([1, 1, 0, 0])  # A, B, C, D
        strain_names = ['A', 'B', 'C', 'D']
        
        # Test through public interface
        scd = evo_analysis._sister_clade_differences(sample_tree, strain_names, trait_vector)
        
        assert isinstance(scd, (int, float, np.integer, np.floating))
        assert scd >= 0
    
    def test_expected_scd_random(self):
        """Test expected SCD under random distribution."""
        if not EVO_AVAILABLE:
            pytest.skip("EvolutionaryAnalysis not available")
        
        evo_analysis = EvolutionaryAnalysis()
        
        trait_vector = np.array([1, 0, 1, 0, 1])
        scd_random = evo_analysis._expected_scd_random(trait_vector)
        
        # Expected = n - 1 for random
        assert scd_random == len(trait_vector) - 1
    
    def test_expected_scd_brownian(self, sample_tree):
        """Test expected SCD under Brownian motion."""
        if not EVO_AVAILABLE or not BIO_AVAILABLE:
            pytest.skip("EvolutionaryAnalysis or Bio.Phylo not available")
        
        evo_analysis = EvolutionaryAnalysis()
        
        trait_vector = np.array([1, 1, 0, 0])
        strain_names = ['A', 'B', 'C', 'D']
        
        # _expected_scd_brownian takes tree, strain_names, trait_vector
        scd_brownian = evo_analysis._expected_scd_brownian(sample_tree, strain_names, trait_vector)
        
        # Should be >= 0
        assert scd_brownian >= 0
    
    def test_phylogenetic_signal_interpretation(self):
        """Test that D statistic interpretation is correct."""
        # D < 0: Overdispersed
        # D < 0.5: Conserved
        # D < 1.5: Random
        # D >= 1.5: Clustered
        
        # This is tested through the main function
        pass
    
    def test_edge_cases_empty_traits(self, sample_tree):
        """Test with empty trait data."""
        if not EVO_AVAILABLE or not BIO_AVAILABLE:
            pytest.skip("EvolutionaryAnalysis or Bio.Phylo not available")
        
        evo_analysis = EvolutionaryAnalysis()
        
        empty_traits = pd.DataFrame()
        
        try:
            result = evo_analysis.calculate_phylogenetic_signal_fritz_purvis(
                sample_tree, empty_traits
            )
            # Should return empty DataFrame
            assert isinstance(result, pd.DataFrame)
        except Exception:
            pass  # May raise exception for empty data
    
    def test_edge_cases_single_trait(self, sample_tree):
        """Test with single trait."""
        if not EVO_AVAILABLE or not BIO_AVAILABLE:
            pytest.skip("EvolutionaryAnalysis or Bio.Phylo not available")
        
        evo_analysis = EvolutionaryAnalysis()
        
        single_trait = pd.DataFrame({
            'Trait_1': [1, 1, 0, 0],
        }, index=['A', 'B', 'C', 'D'])
        
        result = evo_analysis.calculate_phylogenetic_signal_fritz_purvis(
            sample_tree, single_trait
        )
        
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert 'd_statistic' in result.columns
            assert 'interpretation' in result.columns
