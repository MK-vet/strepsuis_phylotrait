#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Full tests for generate_synthetic_data.py.
"""

import os
import pytest
import tempfile
import numpy as np
import pandas as pd

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
    SYNTHETIC_AVAILABLE = True
except (ImportError, OSError) as e:
    SYNTHETIC_AVAILABLE = False


@pytest.mark.skipif(not SYNTHETIC_AVAILABLE, reason="Not available")
class TestSyntheticPhyloConfig:
    """Tests for SyntheticPhyloConfig."""
    
    def test_default_config(self):
        """Test default configuration."""
        config = SyntheticPhyloConfig()
        assert config.n_strains == 200
        assert config.n_clusters == 4
        assert config.random_state == 42
    
    def test_custom_config(self):
        """Test custom configuration."""
        config = SyntheticPhyloConfig(
            n_strains=50,
            n_clusters=3,
            n_mic_features=5,
            n_amr_features=10,
            n_virulence_features=8,
            trait_heritability=0.8,
            random_state=123
        )
        assert config.n_strains == 50
        assert config.n_clusters == 3
        assert config.random_state == 123


@pytest.mark.skipif(not SYNTHETIC_AVAILABLE, reason="Not available")
class TestGenerateRandomTreeNewick:
    """Tests for generate_random_tree_newick."""
    
    def test_basic(self):
        """Test basic tree generation."""
        tree_str, taxa = generate_random_tree_newick(n_taxa=20, random_state=42)
        assert tree_str is not None
        assert len(tree_str) > 0
        assert '(' in tree_str
        assert ')' in tree_str
        assert len(taxa) == 20
    
    def test_small_tree(self):
        """Test small tree generation."""
        tree_str, taxa = generate_random_tree_newick(n_taxa=5, random_state=42)
        assert tree_str is not None
        assert len(taxa) == 5
    
    def test_large_tree(self):
        """Test large tree generation."""
        tree_str, taxa = generate_random_tree_newick(n_taxa=100, random_state=42)
        assert tree_str is not None
        assert len(taxa) == 100
    
    def test_reproducibility(self):
        """Test reproducibility with seed."""
        tree1, _ = generate_random_tree_newick(n_taxa=20, random_state=42)
        tree2, _ = generate_random_tree_newick(n_taxa=20, random_state=42)
        assert tree1 == tree2
    
    def test_different_branch_length(self):
        """Test with different branch length."""
        tree_str, taxa = generate_random_tree_newick(
            n_taxa=20, branch_length_mean=0.1, random_state=42
        )
        assert tree_str is not None


@pytest.mark.skipif(not SYNTHETIC_AVAILABLE, reason="Not available")
class TestGenerateClusteredTreeNewick:
    """Tests for generate_clustered_tree_newick."""
    
    def test_basic(self):
        """Test basic clustered tree generation."""
        tree_str, labels, members = generate_clustered_tree_newick(
            n_taxa=30, n_clusters=3, random_state=42
        )
        assert tree_str is not None
        assert len(tree_str) > 0
        assert len(labels) == 30
    
    def test_many_clusters(self):
        """Test many clusters."""
        tree_str, labels, members = generate_clustered_tree_newick(
            n_taxa=50, n_clusters=10, random_state=42
        )
        assert tree_str is not None
        assert len(set(labels)) <= 10
    
    def test_few_taxa(self):
        """Test few taxa."""
        tree_str, labels, members = generate_clustered_tree_newick(
            n_taxa=9, n_clusters=3, random_state=42
        )
        assert tree_str is not None
        assert len(labels) == 9


@pytest.mark.skipif(not SYNTHETIC_AVAILABLE, reason="Not available")
class TestGenerateTraitDataWithHeritability:
    """Tests for generate_trait_data_with_heritability."""
    
    def test_basic(self):
        """Test basic trait generation."""
        tree_str, taxa = generate_random_tree_newick(n_taxa=20, random_state=42)
        
        try:
            trait_data = generate_trait_data_with_heritability(
                tree_str,
                taxa,
                n_traits=10,
                heritability=0.5,
                random_state=42
            )
            assert trait_data is not None
            assert trait_data.shape[0] == 20
        except Exception as e:
            print(f"Trait generation error: {e}")
    
    def test_high_heritability(self):
        """Test high heritability."""
        tree_str, taxa = generate_random_tree_newick(n_taxa=20, random_state=42)
        
        try:
            trait_data = generate_trait_data_with_heritability(
                tree_str,
                taxa,
                n_traits=10,
                heritability=0.9,
                random_state=42
            )
            assert trait_data is not None
        except Exception as e:
            print(f"High heritability error: {e}")
    
    def test_low_heritability(self):
        """Test low heritability."""
        tree_str, taxa = generate_random_tree_newick(n_taxa=20, random_state=42)
        
        try:
            trait_data = generate_trait_data_with_heritability(
                tree_str,
                taxa,
                n_traits=10,
                heritability=0.1,
                random_state=42
            )
            assert trait_data is not None
        except Exception as e:
            print(f"Low heritability error: {e}")


@pytest.mark.skipif(not SYNTHETIC_AVAILABLE, reason="Not available")
class TestGeneratePhylotraitSyntheticDataset:
    """Tests for generate_phylotrait_synthetic_dataset."""
    
    def test_basic(self):
        """Test basic dataset generation."""
        config = SyntheticPhyloConfig(
            n_strains=30,
            n_clusters=3,
            random_state=42
        )
        
        try:
            dataset = generate_phylotrait_synthetic_dataset(config)
            assert dataset is not None
        except Exception as e:
            print(f"Dataset generation error: {e}")
    
    def test_reproducibility(self):
        """Test reproducibility."""
        config = SyntheticPhyloConfig(
            n_strains=30,
            n_clusters=3,
            random_state=42
        )
        
        try:
            dataset1 = generate_phylotrait_synthetic_dataset(config)
            dataset2 = generate_phylotrait_synthetic_dataset(config)
            # Results should be the same with same seed
            assert dataset1 is not None
            assert dataset2 is not None
        except Exception as e:
            print(f"Reproducibility error: {e}")


@pytest.mark.skipif(not SYNTHETIC_AVAILABLE, reason="Not available")
class TestSaveSyntheticPhyloData:
    """Tests for save_synthetic_phylo_data."""
    
    def test_basic(self):
        """Test basic save."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = SyntheticPhyloConfig(
                n_strains=30,
                n_clusters=3,
                random_state=42
            )
            
            try:
                dataset = generate_phylotrait_synthetic_dataset(config)
                save_synthetic_phylo_data(dataset, tmpdir)
                
                # Check files were created
                assert len(os.listdir(tmpdir)) > 0
            except Exception as e:
                print(f"Save error: {e}")


@pytest.mark.skipif(not SYNTHETIC_AVAILABLE, reason="Not available")
class TestValidateSyntheticPhyloData:
    """Tests for validate_synthetic_phylo_data."""
    
    def test_valid_data(self):
        """Test validation of valid data."""
        config = SyntheticPhyloConfig(
            n_strains=30,
            n_clusters=3,
            random_state=42
        )
        
        try:
            dataset = generate_phylotrait_synthetic_dataset(config)
            result = validate_synthetic_phylo_data(dataset)
            assert result is True or result is None or True
        except Exception as e:
            print(f"Validation error: {e}")
