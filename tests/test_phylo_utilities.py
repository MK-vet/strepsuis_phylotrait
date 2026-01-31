#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for phylogenetic utility functions.

Tests tree loading, trait data loading, and basic phylogenetic calculations.
"""

import os
import tempfile
import pytest
import pandas as pd
import numpy as np

try:
    from Bio import Phylo
    from io import StringIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        load_tree,
        load_trait_data,
    )
    FUNCTIONS_AVAILABLE = True
except ImportError:
    FUNCTIONS_AVAILABLE = False


@pytest.mark.skipif(not FUNCTIONS_AVAILABLE, reason="Functions not available")
class TestLoadTree:
    """Test tree loading."""
    
    @pytest.fixture
    def sample_tree_file(self):
        """Create sample tree file."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        fd, path = tempfile.mkstemp(suffix='.newick')
        with os.fdopen(fd, 'w') as f:
            f.write("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);")
        yield path
        os.remove(path)
    
    def test_load_tree_from_file(self, sample_tree_file):
        """Test loading tree from file."""
        if not BIO_AVAILABLE:
            pytest.skip("Bio.Phylo not available")
        
        tree = load_tree(sample_tree_file)
        
        assert tree is not None
        assert len(list(tree.get_terminals())) == 4


@pytest.mark.skipif(not FUNCTIONS_AVAILABLE, reason="Functions not available")
class TestLoadTraitData:
    """Test trait data loading."""
    
    @pytest.fixture
    def sample_trait_file(self):
        """Create sample trait CSV file."""
        fd, path = tempfile.mkstemp(suffix='.csv')
        df = pd.DataFrame({
            'Strain_ID': ['A', 'B', 'C', 'D'],
            'Trait_1': [1, 1, 0, 0],
            'Trait_2': [0, 0, 1, 1],
        })
        df.to_csv(path, index=False)
        yield path
        os.remove(path)
    
    def test_load_trait_data_from_file(self, sample_trait_file):
        """Test loading trait data from file."""
        trait_data = load_trait_data(sample_trait_file)
        
        assert isinstance(trait_data, pd.DataFrame)
        assert 'Strain_ID' in trait_data.columns or 'Trait_1' in trait_data.columns
