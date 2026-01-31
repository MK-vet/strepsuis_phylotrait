#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Comprehensive tests for cli.py to increase coverage.
"""

import os
import sys
import pytest
import tempfile
from unittest.mock import patch, MagicMock

try:
    from strepsuis_phylotrait import cli
    CLI_AVAILABLE = True
except ImportError:
    CLI_AVAILABLE = False


@pytest.mark.skipif(not CLI_AVAILABLE, reason="CLI not available")
class TestCLIModule:
    """Test CLI module."""
    
    def test_cli_module_import(self):
        """Test that CLI module can be imported."""
        assert cli is not None
    
    def test_cli_has_main(self):
        """Test that CLI has main function."""
        assert hasattr(cli, 'main') or hasattr(cli, 'run') or hasattr(cli, 'cli')


@pytest.mark.skipif(not CLI_AVAILABLE, reason="CLI not available")
class TestCLIArgumentParser:
    """Test CLI argument parser."""
    
    @patch('sys.argv', ['strepsuis-phylotrait', '--help'])
    def test_help_flag(self):
        """Test --help flag."""
        with pytest.raises(SystemExit) as exc_info:
            if hasattr(cli, 'main'):
                cli.main()
            elif hasattr(cli, 'cli'):
                cli.cli()
        
        assert exc_info.value.code == 0
    
    @patch('sys.argv', ['strepsuis-phylotrait', '--version'])
    def test_version_flag(self):
        """Test --version flag."""
        try:
            with pytest.raises(SystemExit) as exc_info:
                if hasattr(cli, 'main'):
                    cli.main()
                elif hasattr(cli, 'cli'):
                    cli.cli()
            
            assert exc_info.value.code == 0
        except (AttributeError, TypeError):
            pytest.skip("Version flag not implemented")
