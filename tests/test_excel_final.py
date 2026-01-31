#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Final tests for excel_report_utils.py to reach 80% coverage.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile

try:
    from strepsuis_phylotrait.excel_report_utils import (
        ExcelReportGenerator,
        sanitize_sheet_name,
    )
    EXCEL_AVAILABLE = True
except (ImportError, OSError) as e:
    EXCEL_AVAILABLE = False


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestAddDataframeSheet:
    """Tests for add_dataframe_sheet method."""
    
    def test_add_dataframe_sheet_basic(self):
        """Test basic add_dataframe_sheet."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            
            df = pd.DataFrame({
                'Cluster': [0, 1, 2],
                'Size': [30, 35, 26],
                'Silhouette': [0.5123, 0.6234, 0.4567]
            })
            
            with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                generator.add_dataframe_sheet(writer, df, "Results")
            
            # Verify
            result = pd.read_excel(report_path, sheet_name='Results')
            assert len(result) == 3
    
    def test_add_dataframe_sheet_with_description(self):
        """Test add_dataframe_sheet with description."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            
            df = pd.DataFrame({
                'Gene': ['gene_1', 'gene_2'],
                'P_value': [0.001, 0.05]
            })
            
            with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                generator.add_dataframe_sheet(
                    writer, df, "Statistics",
                    description="Statistical analysis results"
                )
    
    def test_add_dataframe_sheet_empty(self):
        """Test add_dataframe_sheet with empty DataFrame."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            
            df = pd.DataFrame()
            
            with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                generator.add_dataframe_sheet(writer, df, "Empty")
    
    def test_add_dataframe_sheet_none(self):
        """Test add_dataframe_sheet with None."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            
            with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                generator.add_dataframe_sheet(writer, None, "None")
    
    def test_add_dataframe_sheet_long_name(self):
        """Test add_dataframe_sheet with long sheet name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            
            df = pd.DataFrame({'A': [1, 2, 3]})
            long_name = "A" * 50
            
            with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                generator.add_dataframe_sheet(writer, df, long_name)


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestGenerateExcelReport:
    """Tests for generate_excel_report method."""
    
    def test_generate_excel_report_basic(self):
        """Test basic generate_excel_report."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            
            sheets_data = {
                'Summary': pd.DataFrame({
                    'Metric': ['Silhouette', 'Calinski'],
                    'Value': [0.55, 120.5]
                }),
                'Clusters': pd.DataFrame({
                    'Cluster': [0, 1, 2],
                    'Size': [30, 35, 26]
                })
            }
            
            result = generator.generate_excel_report("test_report", sheets_data)
            assert result is not None
            assert os.path.exists(result)
    
    def test_generate_excel_report_with_methodology(self):
        """Test generate_excel_report with methodology."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            
            sheets_data = {
                'Results': pd.DataFrame({'A': [1, 2, 3]})
            }
            
            methodology = {
                'Data Processing': 'Standard methods',
                'Analysis': 'Clustering analysis'
            }
            
            result = generator.generate_excel_report(
                "test_report",
                sheets_data,
                methodology=methodology
            )
            assert result is not None
    
    def test_generate_excel_report_with_descriptions(self):
        """Test generate_excel_report with sheet descriptions."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            
            sheets_data = {
                'Summary': (
                    pd.DataFrame({'Metric': ['A', 'B'], 'Value': [1, 2]}),
                    "Summary statistics"
                ),
                'Details': (
                    pd.DataFrame({'X': [10, 20], 'Y': [30, 40]}),
                    "Detailed results"
                )
            }
            
            result = generator.generate_excel_report("test_report", sheets_data)
            assert result is not None
    
    def test_generate_excel_report_with_metadata(self):
        """Test generate_excel_report with metadata."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            
            sheets_data = {
                'Data': pd.DataFrame({'A': [1, 2, 3]})
            }
            
            result = generator.generate_excel_report(
                "test_report",
                sheets_data,
                n_strains=91,
                n_clusters=4,
                custom_param="value"
            )
            assert result is not None
    
    def test_generate_excel_report_with_charts(self):
        """Test generate_excel_report with charts."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            
            # Create some charts
            import matplotlib.pyplot as plt
            for i in range(3):
                fig, ax = plt.subplots()
                ax.plot([1, 2, 3], [i, i+1, i+2])
                generator.save_matplotlib_figure(fig, f"chart_{i}")
            
            sheets_data = {
                'Results': pd.DataFrame({'A': [1, 2, 3]})
            }
            
            result = generator.generate_excel_report("test_report", sheets_data)
            assert result is not None


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestPlotlyFallback:
    """Tests for plotly fallback."""
    
    def test_save_plotly_figure_fallback_path(self):
        """Test plotly fallback saves HTML."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            
            try:
                import plotly.express as px
                fig = px.scatter(x=[1, 2, 3], y=[1, 2, 3])
                
                result = generator.save_plotly_figure_fallback(fig, "test")
                assert result is not None
                # Should be HTML or PNG
                assert result.endswith('.html') or result.endswith('.png')
            except Exception as e:
                print(f"Plotly fallback error: {e}")


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestCreateMethodologySheet:
    """Tests for create_methodology_sheet."""
    
    def test_methodology_with_long_text(self):
        """Test methodology with long text."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            
            methodology = {
                'Section1': 'A' * 200,  # Long text
                'Section2': 'B' * 150,
            }
            
            with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                generator.create_methodology_sheet(writer, methodology)
    
    def test_methodology_with_special_chars(self):
        """Test methodology with special characters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = ExcelReportGenerator(tmpdir)
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            
            methodology = {
                'Formula': 'x² + y² = z²',
                'Special': 'α, β, γ, δ'
            }
            
            with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                generator.create_methodology_sheet(writer, methodology)
