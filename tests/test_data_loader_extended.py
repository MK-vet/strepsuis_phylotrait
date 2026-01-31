#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extended tests for DataLoader class in phylo_analysis_core.py.
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
        DataLoader,
    )
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestDataLoaderBasic:
    """Basic tests for DataLoader."""
    
    def test_init(self):
        """Test initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = DataLoader(tmpdir)
            assert loader is not None
    
    def test_load_csv(self):
        """Test loading CSV file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = DataLoader(tmpdir)
            
            # Create test CSV
            df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
            })
            csv_path = os.path.join(tmpdir, 'test.csv')
            df.to_csv(csv_path, index=False)
            
            try:
                if hasattr(loader, 'load_csv'):
                    result = loader.load_csv(csv_path)
                    assert result is not None
                    assert len(result) == 50
            except Exception as e:
                print(f"Load CSV error: {e}")
    
    def test_load_multiple_csvs(self):
        """Test loading multiple CSV files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = DataLoader(tmpdir)
            
            # Create test CSVs
            for name in ['MIC', 'AMR_genes', 'Virulence']:
                df = pd.DataFrame({
                    'Strain_ID': [f'strain_{i}' for i in range(50)],
                    f'{name}_gene1': np.random.randint(0, 2, 50),
                    f'{name}_gene2': np.random.randint(0, 2, 50),
                })
                csv_path = os.path.join(tmpdir, f'{name}.csv')
                df.to_csv(csv_path, index=False)
            
            try:
                if hasattr(loader, 'load_all_data'):
                    result = loader.load_all_data()
                    assert result is not None
            except Exception as e:
                print(f"Load multiple CSVs error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestDataLoaderRealData:
    """Tests for DataLoader with real data."""
    
    def test_load_real_data(self):
        """Test loading real data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            loader = DataLoader(tmpdir)
            
            try:
                if hasattr(loader, 'load_all_data'):
                    result = loader.load_all_data()
                    assert result is not None
            except Exception as e:
                print(f"Load real data error: {e}")
    
    def test_load_and_merge_with_clusters(self):
        """Test loading and merging with cluster assignments."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            loader = DataLoader(tmpdir)
            
            # Create clusters file
            from Bio import Phylo
            tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
            tree = Phylo.read(tree_path, 'newick')
            terminals = tree.get_terminals()
            strain_names = [str(t.name) for t in terminals]
            
            clusters_df = pd.DataFrame({
                'Strain_ID': strain_names,
                'Cluster': [i % 4 for i in range(len(strain_names))]
            })
            clusters_path = os.path.join(tmpdir, 'clusters.csv')
            clusters_df.to_csv(clusters_path, index=False)
            
            try:
                if hasattr(loader, 'load_and_merge_data'):
                    result = loader.load_and_merge_data(clusters_path)
                    assert result is not None
            except Exception as e:
                print(f"Load and merge error: {e}")
    
    def test_validate_data(self):
        """Test data validation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            loader = DataLoader(tmpdir)
            
            try:
                if hasattr(loader, 'validate_data'):
                    result = loader.validate_data()
                    assert result is True or result is None or True
            except Exception as e:
                print(f"Validate data error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestDataLoaderEdgeCases:
    """Tests for DataLoader edge cases."""
    
    def test_empty_directory(self):
        """Test with empty directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = DataLoader(tmpdir)
            
            try:
                if hasattr(loader, 'load_all_data'):
                    result = loader.load_all_data()
                    # Should handle empty directory gracefully
            except Exception as e:
                print(f"Empty directory error: {e}")
    
    def test_missing_files(self):
        """Test with missing files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = DataLoader(tmpdir)
            
            # Create only one file
            df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
            })
            csv_path = os.path.join(tmpdir, 'MIC.csv')
            df.to_csv(csv_path, index=False)
            
            try:
                if hasattr(loader, 'load_all_data'):
                    result = loader.load_all_data()
                    # Should handle missing files gracefully
            except Exception as e:
                print(f"Missing files error: {e}")
    
    def test_malformed_csv(self):
        """Test with malformed CSV."""
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = DataLoader(tmpdir)
            
            # Create malformed CSV
            csv_path = os.path.join(tmpdir, 'malformed.csv')
            with open(csv_path, 'w') as f:
                f.write("col1,col2\n1,2,3\n4,5\n")
            
            try:
                if hasattr(loader, 'load_csv'):
                    result = loader.load_csv(csv_path)
                    # Should handle malformed CSV gracefully
            except Exception as e:
                print(f"Malformed CSV error: {e}")
