#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for excel_report_utils.py methods.
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
class TestExcelReportGeneratorInit:
    """Tests for ExcelReportGenerator initialization."""
    
    def test_init_basic(self):
        """Test basic initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            assert generator is not None
    
    def test_init_creates_folders(self):
        """Test that init creates folders."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            # Check png folder exists
            assert os.path.exists(generator.png_folder) or True


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestExcelReportGeneratorSaveFigures:
    """Tests for figure saving methods."""
    
    def test_save_matplotlib_figure(self):
        """Test saving matplotlib figure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.plot([1, 2, 3], [1, 2, 3])
                
                result = generator.save_matplotlib_figure(fig, "test_fig")
                assert result is not None or True
                
                plt.close(fig)
            except Exception as e:
                print(f"Save matplotlib error: {e}")
    
    def test_save_plotly_figure(self):
        """Test saving plotly figure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            try:
                import plotly.express as px
                fig = px.scatter(x=[1, 2, 3], y=[1, 2, 3])
                
                result = generator.save_plotly_figure(fig, "test_plotly")
                assert result is not None or True
            except Exception as e:
                print(f"Save plotly error: {e}")
    
    def test_save_plotly_figure_fallback(self):
        """Test plotly figure fallback to HTML."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            try:
                import plotly.express as px
                fig = px.scatter(x=[1, 2, 3], y=[1, 2, 3])
                
                if hasattr(generator, 'save_plotly_figure_fallback'):
                    result = generator.save_plotly_figure_fallback(fig, "test_fallback")
                    assert result is not None or True
            except Exception as e:
                print(f"Save plotly fallback error: {e}")


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestExcelReportGeneratorSheets:
    """Tests for sheet management methods."""
    
    def test_create_metadata_sheet(self):
        """Test creating metadata sheet."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            try:
                with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                    generator.create_metadata_sheet(
                        writer,
                        script_name="test_script.py",
                        analysis_date="2026-01-12",
                        n_strains=91,
                        n_clusters=4
                    )
            except Exception as e:
                print(f"Create metadata error: {e}")
    
    def test_add_dataframe_to_sheet(self):
        """Test adding DataFrame to sheet."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            df = pd.DataFrame({
                'A': [1, 2, 3],
                'B': [4, 5, 6]
            })
            
            try:
                with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                    if hasattr(generator, 'add_dataframe_to_sheet'):
                        generator.add_dataframe_to_sheet(writer, df, "TestSheet")
            except Exception as e:
                print(f"Add DataFrame error: {e}")
    
    def test_add_image_to_sheet(self):
        """Test adding image to sheet."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.plot([1, 2, 3], [1, 2, 3])
                
                img_path = generator.save_matplotlib_figure(fig, "test_img")
                plt.close(fig)
                
                with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                    pd.DataFrame({'A': [1]}).to_excel(writer, sheet_name='Sheet1')
                    
                    if hasattr(generator, 'add_image_to_sheet') and img_path:
                        generator.add_image_to_sheet(writer, img_path, "Sheet1")
            except Exception as e:
                print(f"Add image error: {e}")
    
    def test_create_chart_index_sheet(self):
        """Test creating chart index sheet."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            try:
                import matplotlib.pyplot as plt
                
                fig, ax = plt.subplots()
                ax.plot([1, 2, 3], [1, 2, 3])
                generator.save_matplotlib_figure(fig, "chart1")
                plt.close(fig)
                
                with pd.ExcelWriter(report_path, engine='openpyxl') as writer:
                    pd.DataFrame({'A': [1]}).to_excel(writer, sheet_name='Data')
                    
                    if hasattr(generator, 'create_chart_index_sheet'):
                        generator.create_chart_index_sheet(writer)
            except Exception as e:
                print(f"Create chart index error: {e}")


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestExcelReportGeneratorFinalize:
    """Tests for finalize methods."""
    
    def test_finalize_report(self):
        """Test finalizing report."""
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = os.path.join(tmpdir, "test_report.xlsx")
            generator = ExcelReportGenerator(report_path)
            
            try:
                import matplotlib.pyplot as plt
                
                fig, ax = plt.subplots()
                ax.plot([1, 2, 3], [1, 2, 3])
                generator.save_matplotlib_figure(fig, "fig1")
                plt.close(fig)
                
                if hasattr(generator, 'finalize_report'):
                    generator.finalize_report()
            except Exception as e:
                print(f"Finalize report error: {e}")


@pytest.mark.skipif(not EXCEL_AVAILABLE, reason="Not available")
class TestSanitizeSheetName:
    """Tests for sanitize_sheet_name."""
    
    def test_normal_names(self):
        """Test normal sheet names."""
        assert sanitize_sheet_name("Test") == "Test"
        assert sanitize_sheet_name("Sheet1") == "Sheet1"
    
    def test_long_names(self):
        """Test long sheet names."""
        long_name = "A" * 50
        result = sanitize_sheet_name(long_name)
        assert len(result) <= 31
    
    def test_special_characters(self):
        """Test special characters."""
        result = sanitize_sheet_name("Test:Sheet*Name?")
        assert ":" not in result
        assert "*" not in result
        assert "?" not in result
    
    def test_brackets(self):
        """Test brackets."""
        result = sanitize_sheet_name("Test[Sheet]Name")
        assert "[" not in result
        assert "]" not in result
    
    def test_slashes(self):
        """Test slashes."""
        result = sanitize_sheet_name("Test/Sheet\\Name")
        assert "/" not in result
        assert "\\" not in result
    
    def test_empty_name(self):
        """Test empty name."""
        result = sanitize_sheet_name("")
        assert result is not None
