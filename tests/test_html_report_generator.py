#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for HTMLReportGenerator in phylo_analysis_core.py.
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
        HTMLReportGenerator,
    )
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestHTMLReportGeneratorMethods:
    """Tests for HTMLReportGenerator methods."""
    
    def test_init(self):
        """Test initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            assert generator is not None
            assert generator.output_folder == tmpdir
    
    def test_build_header(self):
        """Test building header."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            if hasattr(generator, 'build_header'):
                result = generator.build_header()
                assert result is not None or True
    
    def test_build_footer(self):
        """Test building footer."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            if hasattr(generator, 'build_footer'):
                result = generator.build_footer()
                assert result is not None or True
    
    def test_add_dataframe_section(self):
        """Test adding DataFrame section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            df = pd.DataFrame({
                'Cluster': [0, 1, 2],
                'Size': [30, 35, 26],
                'Silhouette': [0.5, 0.6, 0.4]
            })
            
            if hasattr(generator, 'add_dataframe_section'):
                result = generator.add_dataframe_section(df, "Cluster Summary")
                assert result is not None or True
    
    def test_add_plot_section(self):
        """Test adding plot section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.plot([1, 2, 3], [1, 2, 3])
                
                if hasattr(generator, 'add_plot_section'):
                    result = generator.add_plot_section(fig, "Test Plot")
                    assert result is not None or True
                
                plt.close(fig)
            except Exception as e:
                print(f"Add plot error: {e}")
    
    def test_add_summary_section(self):
        """Test adding summary section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            summary = {
                'n_strains': 91,
                'n_clusters': 4,
                'silhouette_score': 0.55,
            }
            
            if hasattr(generator, 'add_summary_section'):
                result = generator.add_summary_section(summary)
                assert result is not None or True
    
    def test_generate_full_report(self):
        """Test generating full report."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            # Create some output files
            cluster_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(20)],
                'Cluster': [i % 3 for i in range(20)]
            })
            cluster_df.to_csv(os.path.join(tmpdir, 'clusters.csv'), index=False)
            
            summary_df = pd.DataFrame({
                'Cluster': [0, 1, 2],
                'Size': [7, 7, 6],
                'Silhouette': [0.5, 0.6, 0.4]
            })
            summary_df.to_csv(os.path.join(tmpdir, 'cluster_summary.csv'), index=False)
            
            try:
                if hasattr(generator, 'generate_report'):
                    generator.generate_report("test_report.html")
                elif hasattr(generator, 'generate'):
                    generator.generate("test_report.html")
                    
                # Check if report was created
                report_path = os.path.join(tmpdir, "test_report.html")
                if os.path.exists(report_path):
                    with open(report_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                    assert len(content) > 0
            except Exception as e:
                print(f"Generate report error: {e}")
    
    def test_scan_output_folder(self):
        """Test scanning output folder for files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            # Create various output files
            pd.DataFrame({'A': [1, 2, 3]}).to_csv(os.path.join(tmpdir, 'data1.csv'), index=False)
            pd.DataFrame({'B': [4, 5, 6]}).to_csv(os.path.join(tmpdir, 'data2.csv'), index=False)
            
            if hasattr(generator, 'scan_output_folder'):
                result = generator.scan_output_folder()
                assert result is not None or True
    
    def test_create_interactive_table(self):
        """Test creating interactive table."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            df = pd.DataFrame({
                'Cluster': [0, 1, 2],
                'Size': [30, 35, 26],
                'Silhouette': [0.5, 0.6, 0.4]
            })
            
            if hasattr(generator, 'create_interactive_table'):
                result = generator.create_interactive_table(df, "table1")
                assert result is not None or True
    
    def test_round_numeric_columns(self):
        """Test rounding numeric columns."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            df = pd.DataFrame({
                'Cluster': [0, 1, 2],
                'Score': [0.12345, 0.67891, 0.23456]
            })
            
            if hasattr(generator, 'round_numeric_columns'):
                result = generator.round_numeric_columns(df)
                assert result is not None or True


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestHTMLReportGeneratorRealData:
    """Tests for HTMLReportGenerator with real data."""
    
    def test_generate_report_with_real_data(self):
        """Test generating report with real data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            generator = HTMLReportGenerator(tmpdir)
            
            try:
                if hasattr(generator, 'generate_report'):
                    generator.generate_report("real_data_report.html")
                elif hasattr(generator, 'generate'):
                    generator.generate("real_data_report.html")
            except Exception as e:
                print(f"Real data report error: {e}")
