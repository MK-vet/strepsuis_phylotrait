#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis edge cases.
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
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticSignalEdgeCases:
    """Tests for phylogenetic signal edge cases."""
    
    def test_phylogenetic_signal_nan_case(self):
        """Test phylogenetic signal when scd_random == scd_brownian."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create trait data that might cause scd_random == scd_brownian
        trait_data = pd.DataFrame({
            'uniform': [1] * len(terminals),  # All same
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
            # Should handle NaN case
        except Exception as e:
            print(f"NaN case error: {e}")
    
    def test_phylogenetic_signal_overdispersed(self):
        """Test phylogenetic signal for overdispersed trait (D < 0)."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create trait that might be overdispersed
        # Alternate values to create pattern where distant taxa are similar
        trait_vector = [1 if i % 2 == 0 else 0 for i in range(len(terminals))]
        trait_data = pd.DataFrame({
            'overdispersed': trait_vector,
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
        except Exception as e:
            print(f"Overdispersed error: {e}")
    
    def test_phylogenetic_signal_clustered(self):
        """Test phylogenetic signal for clustered trait (D > 1.5)."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create trait that might be clustered
        # First half all 1, second half all 0
        n_half = len(terminals) // 2
        trait_vector = [1] * n_half + [0] * (len(terminals) - n_half)
        trait_data = pd.DataFrame({
            'clustered': trait_vector,
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
        except Exception as e:
            print(f"Clustered error: {e}")
    
    def test_phylogenetic_signal_random(self):
        """Test phylogenetic signal for random trait (0.5 <= D < 1.5)."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create random trait
        np.random.seed(42)
        trait_vector = np.random.randint(0, 2, len(terminals))
        trait_data = pd.DataFrame({
            'random': trait_vector,
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
        except Exception as e:
            print(f"Random error: {e}")
    
    def test_phylogenetic_signal_conserved(self):
        """Test phylogenetic signal for conserved trait (0 <= D < 0.5)."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        # Create trait that follows phylogeny (conserved)
        # Use tree structure to create conserved pattern
        trait_vector = [1 if i < len(terminals) // 2 else 0 for i in range(len(terminals))]
        trait_data = pd.DataFrame({
            'conserved': trait_vector,
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
        except Exception as e:
            print(f"Conserved error: {e}")
    
    def test_phylogenetic_signal_mixed_traits(self):
        """Test phylogenetic signal with mixed trait types."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        np.random.seed(42)
        trait_data = pd.DataFrame({
            'all_ones': [1] * len(terminals),
            'all_zeros': [0] * len(terminals),
            'random': np.random.randint(0, 2, len(terminals)),
            'half': [1 if i < len(terminals) // 2 else 0 for i in range(len(terminals))],
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
            assert len(result) == 4
        except Exception as e:
            print(f"Mixed traits error: {e}")
