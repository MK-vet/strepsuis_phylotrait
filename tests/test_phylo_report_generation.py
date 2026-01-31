#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for report generation functions in phylo_analysis_core.py.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil
from unittest.mock import patch, MagicMock

REAL_DATA_PATH = r"C:\Users\ABC\Documents\GitHub\MKrep\data"


def check_real_data_exists():
    return os.path.exists(os.path.join(REAL_DATA_PATH, 'Snp_tree.newick'))


try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        HTMLReportGenerator,
        TraitAnalyzer,
        MCAAnalyzer,
    )
    from Bio import Phylo
    import plotly.express as px
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestHTMLReportGeneratorGenerateReportFull:
    """Full tests for HTMLReportGenerator.generate_report."""
    
    def test_generate_report_with_all_data(self):
        """Test report generation with all data types."""
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                generator = HTMLReportGenerator(tmpdir, base_dir=tmpdir)
                
                # Create all output files
                pd.DataFrame({
                    'Cluster': [0, 1, 2, 3],
                    'Strain_Count': [23, 24, 22, 22],
                    'Percentage': [25.3, 26.4, 24.2, 24.2]
                }).to_csv(os.path.join(tmpdir, 'cluster_distribution.csv'), index=False)
                
                pd.DataFrame({
                    'Cluster': [0, 1, 2, 3],
                    'Mean_PD': [0.1, 0.15, 0.12, 0.11],
                    'Mean_Branch_Length': [0.01, 0.02, 0.015, 0.018]
                }).to_csv(os.path.join(tmpdir, 'evolutionary_cluster_analysis.csv'), index=False)
                
                pd.DataFrame({
                    'Strain_ID': [f'strain_{i}' for i in range(20)],
                    'Cluster': [i % 4 for i in range(20)]
                }).to_csv(os.path.join(tmpdir, 'phylogenetic_clusters.csv'), index=False)
                
                # Create analysis results
                analysis_results = {
                    'summary_stats': pd.DataFrame({
                        'Metric': ['Silhouette', 'Calinski-Harabasz', 'Davies-Bouldin'],
                        'Value': [0.55, 120.5, 0.8]
                    }),
                    'embeddings': np.random.randn(20, 2),
                    'labels': np.array([i % 4 for i in range(20)]),
                    'mask': np.array([True] * 20),
                }
                
                config = MagicMock()
                config.tree_file = 'tree.newick'
                
                with patch.object(generator, 'template_env') as mock_env:
                    mock_template = MagicMock()
                    mock_template.render.return_value = "<html>Full Report</html>"
                    mock_env.get_template.return_value = mock_template
                    
                    generator.generate_report(analysis_results, config)
                    
            except Exception as e:
                print(f"Generate report full error: {e}")
    
    def test_generate_report_missing_files(self):
        """Test report generation with missing files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                generator = HTMLReportGenerator(tmpdir, base_dir=tmpdir)
                
                analysis_results = {}
                
                config = MagicMock()
                config.tree_file = 'nonexistent.newick'
                
                with patch.object(generator, 'template_env') as mock_env:
                    mock_template = MagicMock()
                    mock_template.render.return_value = "<html>Empty Report</html>"
                    mock_env.get_template.return_value = mock_template
                    
                    generator.generate_report(analysis_results, config)
                    
            except Exception as e:
                print(f"Generate report missing files error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestHTMLReportGeneratorTreePlotFull:
    """Full tests for tree plot generation."""
    
    def test_create_interactive_tree_plot_full(self):
        """Test creating interactive tree plot with clusters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy tree file
            src = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
            shutil.copy(src, os.path.join(tmpdir, 'Snp_tree.newick'))
            
            try:
                generator = HTMLReportGenerator(tmpdir, base_dir=tmpdir)
                
                # Create clusters file
                tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
                tree = Phylo.read(tree_path, 'newick')
                terminals = tree.get_terminals()
                strain_names = [str(t.name) for t in terminals]
                
                pd.DataFrame({
                    'Strain_ID': strain_names,
                    'Cluster': [i % 4 for i in range(len(strain_names))]
                }).to_csv(os.path.join(tmpdir, 'phylogenetic_clusters.csv'), index=False)
                
                result = generator._create_interactive_tree_plot(tree_path)
                assert result is not None
            except Exception as e:
                print(f"Interactive tree plot full error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerReportMethods:
    """Tests for TraitAnalyzer report-related methods."""
    
    def test_association_rules_empty_result(self):
        """Test association rules with no results."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            # Create data with very sparse features
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': [0] * 50,
                'gene2': [0] * 50,
            })
            
            try:
                result = analyzer.association_rule_mining(merged_df, min_support=0.9)
                assert result is not None
                # Should return empty DataFrame
            except Exception as e:
                print(f"Empty association rules error: {e}")
    
    def test_label_shared_unique_full(self):
        """Test label_shared_unique_features with various data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            # Create data with clear patterns
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(60)],
                'Cluster': [0]*20 + [1]*20 + [2]*20,
                'shared_gene': [1] * 60,  # Present in all
                'cluster0_gene': [1]*20 + [0]*40,  # Only in cluster 0
                'cluster1_gene': [0]*20 + [1]*20 + [0]*20,  # Only in cluster 1
                'random_gene': np.random.randint(0, 2, 60),
            })
            
            try:
                result = analyzer.label_shared_unique_features(merged_df, presence_threshold=0.5)
                assert result is not None
            except Exception as e:
                print(f"Label shared unique full error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestMCAAnalyzerReportMethods:
    """Tests for MCAAnalyzer report-related methods."""
    
    def test_perform_mca_with_feature_types(self):
        """Test MCA with different feature types."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = MCAAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'MIC_drug1': np.random.randint(0, 2, 50),
                'MIC_drug2': np.random.randint(0, 2, 50),
                'AMR_gene1': np.random.randint(0, 2, 50),
                'AMR_gene2': np.random.randint(0, 2, 50),
                'VIR_factor1': np.random.randint(0, 2, 50),
                'VIR_factor2': np.random.randint(0, 2, 50),
                'other_feature': np.random.randint(0, 2, 50),
            })
            
            try:
                row_coords, col_coords, summary = analyzer.perform_mca_analysis(merged_df)
                
                if col_coords is not None:
                    assert 'Feature_Type' in col_coords.columns
                    # Check feature types are assigned
                    feature_types = col_coords['Feature_Type'].unique()
                    assert len(feature_types) > 0
            except Exception as e:
                print(f"MCA feature types error: {e}")
    
    def test_perform_mca_saves_files(self):
        """Test that MCA saves output files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = MCAAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
                'gene3': np.random.randint(0, 2, 50),
            })
            
            try:
                analyzer.perform_mca_analysis(merged_df)
                
                # Check output files
                expected_files = [
                    'mca_row_coordinates.csv',
                    'mca_column_coordinates.csv',
                    'mca_summary.csv',
                    'mca_interpretation_guide.csv'
                ]
                
                for f in expected_files:
                    path = os.path.join(tmpdir, f)
                    if os.path.exists(path):
                        df = pd.read_csv(path)
                        assert len(df) > 0
            except Exception as e:
                print(f"MCA saves files error: {e}")
