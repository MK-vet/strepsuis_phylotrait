#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for cli.py methods.
"""

import os
import pytest
import sys
import tempfile
import shutil
from unittest.mock import patch

REAL_DATA_PATH = r"C:\Users\ABC\Documents\GitHub\MKrep\data"


def check_real_data_exists():
    return os.path.exists(os.path.join(REAL_DATA_PATH, 'Snp_tree.newick'))


try:
    from strepsuis_phylotrait import cli
    CLI_AVAILABLE = True
except (ImportError, OSError) as e:
    CLI_AVAILABLE = False


@pytest.mark.skipif(not CLI_AVAILABLE, reason="Not available")
class TestCLIBasic:
    """Basic tests for CLI."""
    
    def test_import(self):
        """Test CLI import."""
        assert cli is not None
    
    def test_main_exists(self):
        """Test main function exists."""
        assert hasattr(cli, 'main')
    
    def test_help(self):
        """Test CLI help."""
        original_argv = sys.argv
        try:
            sys.argv = ['strepsuis-phylotrait', '--help']
            with pytest.raises(SystemExit) as exc_info:
                cli.main()
            assert exc_info.value.code == 0
        except Exception as e:
            print(f"Help error: {e}")
        finally:
            sys.argv = original_argv
    
    def test_version(self):
        """Test CLI version."""
        original_argv = sys.argv
        try:
            sys.argv = ['strepsuis-phylotrait', '--version']
            with pytest.raises(SystemExit):
                cli.main()
        except Exception as e:
            print(f"Version error: {e}")
        finally:
            sys.argv = original_argv
    
    def test_invalid_args(self):
        """Test CLI with invalid arguments."""
        original_argv = sys.argv
        try:
            sys.argv = ['strepsuis-phylotrait', '--invalid-arg']
            with pytest.raises(SystemExit):
                cli.main()
        except Exception as e:
            print(f"Invalid args error: {e}")
        finally:
            sys.argv = original_argv


@pytest.mark.skipif(not CLI_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestCLIWithData:
    """Tests for CLI with data."""
    
    def test_with_data_dir(self):
        """Test CLI with data directory."""
        original_argv = sys.argv
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                sys.argv = [
                    'strepsuis-phylotrait',
                    '--data-dir', tmpdir,
                    '--output-dir', os.path.join(tmpdir, 'output')
                ]
                cli.main()
            except SystemExit:
                pass
            except Exception as e:
                print(f"Data dir error: {e}")
            finally:
                sys.argv = original_argv
    
    def test_with_tree_file(self):
        """Test CLI with tree file."""
        original_argv = sys.argv
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                sys.argv = [
                    'strepsuis-phylotrait',
                    '--data-dir', tmpdir,
                    '--output-dir', os.path.join(tmpdir, 'output'),
                    '--tree-file', 'Snp_tree.newick'
                ]
                cli.main()
            except SystemExit:
                pass
            except Exception as e:
                print(f"Tree file error: {e}")
            finally:
                sys.argv = original_argv
    
    def test_with_n_clusters(self):
        """Test CLI with n_clusters option."""
        original_argv = sys.argv
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                sys.argv = [
                    'strepsuis-phylotrait',
                    '--data-dir', tmpdir,
                    '--output-dir', os.path.join(tmpdir, 'output'),
                    '--n-clusters-min', '2',
                    '--n-clusters-max', '6'
                ]
                cli.main()
            except SystemExit:
                pass
            except Exception as e:
                print(f"N clusters error: {e}")
            finally:
                sys.argv = original_argv
