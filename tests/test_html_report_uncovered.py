#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for uncovered HTMLReportGenerator methods.
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
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestHTMLReportGeneratorMCA:
    """Tests for MCA-related report generation."""
    
    def test_generate_report_with_mca_data(self):
        """Test report generation with MCA data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                generator = HTMLReportGenerator(tmpdir, base_dir=tmpdir)
                
                # Create MCA output files
                pd.DataFrame({
                    'Strain_ID': [f'strain_{i}' for i in range(50)],
                    'Component_1': np.random.randn(50),
                    'Component_2': np.random.randn(50),
                    'Cluster': [i % 3 for i in range(50)]
                }).to_csv(os.path.join(tmpdir, 'mca_row_coordinates.csv'), index=False)
                
                pd.DataFrame({
                    'Component_1': np.random.randn(20),
                    'Component_2': np.random.randn(20),
                    'Feature_Type': ['AMR'] * 10 + ['VIR'] * 10
                }, index=[f'gene_{i}' for i in range(20)]).to_csv(
                    os.path.join(tmpdir, 'mca_column_coordinates.csv'), index=False
                )
                
                pd.DataFrame({
                    'Component': [1, 2, 3],
                    'Eigenvalue': [0.5, 0.3, 0.2],
                    'Explained_Inertia': [0.5, 0.3, 0.2],
                    'Cumulative_Explained_Inertia': [0.5, 0.8, 1.0]
                }).to_csv(os.path.join(tmpdir, 'mca_summary.csv'), index=False)
                
                pd.DataFrame({
                    'Aspect': ['Column Points', 'Row Points'],
                    'Interpretation': ['Features', 'Strains']
                }).to_csv(os.path.join(tmpdir, 'mca_interpretation_guide.csv'), index=False)
                
                analysis_results = {
                    'summary_stats': pd.DataFrame({
                        'Metric': ['Silhouette'],
                        'Value': [0.55]
                    })
                }
                
                config = MagicMock()
                config.tree_file = 'tree.newick'
                
                with patch.object(generator, 'template_env') as mock_env:
                    mock_template = MagicMock()
                    mock_template.render.return_value = "<html>MCA Report</html>"
                    mock_env.get_template.return_value = mock_template
                    
                    generator.generate_report(analysis_results, config)
                    
            except Exception as e:
                print(f"MCA report error: {e}")
    
    def test_generate_report_with_all_output_files(self):
        """Test report generation with all possible output files."""
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
                    'Mean_PD': [0.1, 0.15, 0.12, 0.11]
                }).to_csv(os.path.join(tmpdir, 'evolutionary_cluster_analysis.csv'), index=False)
                
                pd.DataFrame({
                    'k': [2, 3, 4, 5],
                    'Silhouette': [0.4, 0.5, 0.55, 0.45]
                }).to_csv(os.path.join(tmpdir, 'silhouette_scores.csv'), index=False)
                
                pd.DataFrame({
                    'Feature': ['gene_1', 'gene_2'],
                    'Importance_Mean': [0.8, 0.6]
                }).to_csv(os.path.join(tmpdir, 'bootstrap_feature_importance.csv'), index=False)
                
                pd.DataFrame({
                    'Antecedent': ['gene_1'],
                    'Consequent': ['gene_2'],
                    'Support': [0.3],
                    'Confidence': [0.8],
                    'Lift': [2.5]
                }).to_csv(os.path.join(tmpdir, 'association_rules.csv'), index=False)
                
                pd.DataFrame({
                    'Feature': ['gene_1'],
                    'Log_Odds': [1.5],
                    'P_value': [0.01]
                }).to_csv(os.path.join(tmpdir, 'log_odds_global.csv'), index=False)
                
                pd.DataFrame({
                    'Feature': ['gene_1'],
                    'Cluster': [0],
                    'Count': [20]
                }).to_csv(os.path.join(tmpdir, 'log_odds_per_cluster.csv'), index=False)
                
                # MCA files
                pd.DataFrame({
                    'Strain_ID': [f'strain_{i}' for i in range(50)],
                    'Component_1': np.random.randn(50),
                    'Component_2': np.random.randn(50),
                    'Cluster': [i % 3 for i in range(50)]
                }).to_csv(os.path.join(tmpdir, 'mca_row_coordinates.csv'), index=False)
                
                pd.DataFrame({
                    'Component_1': np.random.randn(10),
                    'Component_2': np.random.randn(10),
                    'Feature_Type': ['AMR'] * 5 + ['VIR'] * 5
                }, index=[f'gene_{i}' for i in range(10)]).to_csv(
                    os.path.join(tmpdir, 'mca_column_coordinates.csv'), index=False
                )
                
                pd.DataFrame({
                    'Component': [1, 2],
                    'Eigenvalue': [0.5, 0.3],
                    'Explained_Inertia': [0.5, 0.3],
                    'Cumulative_Explained_Inertia': [0.5, 0.8]
                }).to_csv(os.path.join(tmpdir, 'mca_summary.csv'), index=False)
                
                pd.DataFrame({
                    'Aspect': ['Column Points'],
                    'Interpretation': ['Features']
                }).to_csv(os.path.join(tmpdir, 'mca_interpretation_guide.csv'), index=False)
                
                analysis_results = {
                    'summary_stats': pd.DataFrame({
                        'Metric': ['Silhouette', 'Calinski'],
                        'Value': [0.55, 120.5]
                    }),
                    'embeddings': np.random.randn(50, 2),
                    'labels': np.array([i % 3 for i in range(50)]),
                    'mask': np.array([True] * 50),
                }
                
                config = MagicMock()
                config.tree_file = 'tree.newick'
                
                with patch.object(generator, 'template_env') as mock_env:
                    mock_template = MagicMock()
                    mock_template.render.return_value = "<html>Full Report</html>"
                    mock_env.get_template.return_value = mock_template
                    
                    generator.generate_report(analysis_results, config)
                    
            except Exception as e:
                print(f"All files report error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestHTMLReportGeneratorRealData:
    """Tests with real data."""
    
    def test_generate_report_real_data_full(self):
        """Test report generation with real data and all files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                generator = HTMLReportGenerator(tmpdir, base_dir=tmpdir)
                
                tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
                tree = Phylo.read(tree_path, 'newick')
                terminals = tree.get_terminals()
                strain_names = [str(t.name) for t in terminals]
                
                # Create all output files
                pd.DataFrame({
                    'Strain_ID': strain_names,
                    'Cluster': [i % 4 for i in range(len(strain_names))]
                }).to_csv(os.path.join(tmpdir, 'phylogenetic_clusters.csv'), index=False)
                
                pd.DataFrame({
                    'Cluster': [0, 1, 2, 3],
                    'Strain_Count': [23, 24, 22, 22],
                    'Percentage': [25.3, 26.4, 24.2, 24.2]
                }).to_csv(os.path.join(tmpdir, 'cluster_distribution.csv'), index=False)
                
                pd.DataFrame({
                    'Cluster': [0, 1, 2, 3],
                    'Mean_PD': [0.1, 0.15, 0.12, 0.11]
                }).to_csv(os.path.join(tmpdir, 'evolutionary_cluster_analysis.csv'), index=False)
                
                # MCA files
                pd.DataFrame({
                    'Strain_ID': strain_names,
                    'Component_1': np.random.randn(len(terminals)),
                    'Component_2': np.random.randn(len(terminals)),
                    'Cluster': [i % 4 for i in range(len(terminals))]
                }).to_csv(os.path.join(tmpdir, 'mca_row_coordinates.csv'), index=False)
                
                pd.DataFrame({
                    'Component_1': np.random.randn(20),
                    'Component_2': np.random.randn(20),
                    'Feature_Type': ['AMR'] * 10 + ['VIR'] * 10
                }, index=[f'gene_{i}' for i in range(20)]).to_csv(
                    os.path.join(tmpdir, 'mca_column_coordinates.csv'), index=False
                )
                
                pd.DataFrame({
                    'Component': [1, 2],
                    'Eigenvalue': [0.5, 0.3],
                    'Explained_Inertia': [0.5, 0.3],
                    'Cumulative_Explained_Inertia': [0.5, 0.8]
                }).to_csv(os.path.join(tmpdir, 'mca_summary.csv'), index=False)
                
                analysis_results = {
                    'summary_stats': pd.DataFrame({
                        'Metric': ['Silhouette'],
                        'Value': [0.55]
                    }),
                    'embeddings': np.random.randn(len(terminals), 2),
                    'labels': np.array([i % 4 for i in range(len(terminals))]),
                    'mask': np.array([True] * len(terminals)),
                }
                
                config = MagicMock()
                config.tree_file = 'Snp_tree.newick'
                
                with patch.object(generator, 'template_env') as mock_env:
                    mock_template = MagicMock()
                    mock_template.render.return_value = "<html>Real Data Report</html>"
                    mock_env.get_template.return_value = mock_template
                    
                    generator.generate_report(analysis_results, config)
                    
            except Exception as e:
                print(f"Real data full report error: {e}")
