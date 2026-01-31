#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for generate_synthetic_data.py to increase coverage.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile

try:
    from strepsuis_phylotrait.generate_synthetic_data import (
        SyntheticPhyloConfig,
        SyntheticPhyloMetadata,
        generate_random_tree_newick,
        generate_clustered_tree_newick,
        generate_trait_data_with_heritability,
        generate_phylotrait_synthetic_dataset,
        save_synthetic_phylo_data,
        validate_synthetic_phylo_data,
    )
    SYNTH_AVAILABLE = True
except ImportError:
    SYNTH_AVAILABLE = False


@pytest.mark.skipif(not SYNTH_AVAILABLE, reason="Synthetic data module not available")
class TestSyntheticPhyloConfig:
    """Test SyntheticPhyloConfig class."""
    
    def test_config_default(self):
        """Test default configuration."""
        config = SyntheticPhyloConfig()
        assert config is not None
        assert config.n_taxa > 0
    
    def test_config_custom(self):
        """Test custom configuration."""
        config = SyntheticPhyloConfig(n_taxa=50, n_traits=20)
        assert config.n_taxa == 50
        assert config.n_traits == 20


@pytest.mark.skipif(not SYNTH_AVAILABLE, reason="Synthetic data module not available")
class TestGenerateRandomTreeNewick:
    """Test random tree generation."""
    
    def test_generate_tree_basic(self):
        """Test basic tree generation."""
        tree_str, _ = generate_random_tree_newick(n_taxa=20, random_state=42)
        
        assert tree_str is not None
        assert isinstance(tree_str, str)
        assert len(tree_str) > 0
    
    def test_generate_tree_reproducible(self):
        """Test reproducibility."""
        tree1, _ = generate_random_tree_newick(n_taxa=10, random_state=42)
        tree2, _ = generate_random_tree_newick(n_taxa=10, random_state=42)
        
        assert tree1 == tree2


@pytest.mark.skipif(not SYNTH_AVAILABLE, reason="Synthetic data module not available")
class TestGenerateClusteredTreeNewick:
    """Test clustered tree generation."""
    
    def test_generate_clustered_tree_basic(self):
        """Test basic clustered tree generation."""
        tree_str, _, _ = generate_clustered_tree_newick(
            n_taxa=20, n_clusters=3, random_state=42
        )
        
        assert tree_str is not None
        assert isinstance(tree_str, str)


@pytest.mark.skipif(not SYNTH_AVAILABLE, reason="Synthetic data module not available")
class TestGenerateTraitDataWithHeritability:
    """Test trait data generation."""
    
    def test_generate_traits_basic(self):
        """Test basic trait generation."""
        tree_str, _ = generate_random_tree_newick(n_taxa=20, random_state=42)
        
        traits = generate_trait_data_with_heritability(
            tree_newick=tree_str,
            n_traits=10,
            random_state=42
        )
        
        assert isinstance(traits, pd.DataFrame)
        assert len(traits) == 20
    
    def test_generate_traits_binary(self):
        """Test that traits are binary."""
        tree_str, _ = generate_random_tree_newick(n_taxa=15, random_state=42)
        
        traits = generate_trait_data_with_heritability(
            tree_newick=tree_str,
            n_traits=5,
            random_state=42
        )
        
        for col in traits.columns:
            if col != 'Strain_ID' and col != 'taxon':
                assert set(traits[col].unique()).issubset({0, 1})


@pytest.mark.skipif(not SYNTH_AVAILABLE, reason="Synthetic data module not available")
class TestGeneratePhylotraitSyntheticDataset:
    """Test full synthetic dataset generation."""
    
    def test_generate_dataset_basic(self):
        """Test basic dataset generation."""
        config = SyntheticPhyloConfig(n_taxa=20, n_traits=10)
        
        data, metadata = generate_phylotrait_synthetic_dataset(config)
        
        assert data is not None
        assert metadata is not None
    
    def test_generate_dataset_reproducible(self):
        """Test reproducibility."""
        config1 = SyntheticPhyloConfig(n_taxa=15, n_traits=5, random_state=42)
        config2 = SyntheticPhyloConfig(n_taxa=15, n_traits=5, random_state=42)
        
        data1, _ = generate_phylotrait_synthetic_dataset(config1)
        data2, _ = generate_phylotrait_synthetic_dataset(config2)
        
        pd.testing.assert_frame_equal(data1, data2)


@pytest.mark.skipif(not SYNTH_AVAILABLE, reason="Synthetic data module not available")
class TestSaveSyntheticPhyloData:
    """Test saving synthetic data."""
    
    def test_save_data_basic(self):
        """Test basic data saving."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = SyntheticPhyloConfig(n_taxa=15, n_traits=5)
            data, metadata = generate_phylotrait_synthetic_dataset(config)
            
            result = save_synthetic_phylo_data(data, metadata, tmpdir)
            
            assert isinstance(result, dict)
            files = os.listdir(tmpdir)
            assert len(files) > 0


@pytest.mark.skipif(not SYNTH_AVAILABLE, reason="Synthetic data module not available")
class TestValidateSyntheticPhyloData:
    """Test validation of synthetic data."""
    
    def test_validate_valid_data(self):
        """Test validation of valid data."""
        config = SyntheticPhyloConfig(n_taxa=15, n_traits=5)
        data, metadata = generate_phylotrait_synthetic_dataset(config)
        
        result = validate_synthetic_phylo_data(data, metadata)
        
        assert result is not None
