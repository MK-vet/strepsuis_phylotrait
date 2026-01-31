# Changelog

All notable changes to StrepSuis-PhyloTrait will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2025-11-20

### Added
- Initial release of StrepSuis-PhyloTrait
- Phylogenetic and binary trait analysis
- Tree-aware clustering with phylogenetic constraints
- Evolutionary metrics (phylogenetic diversity, beta diversity)
- Phylogenetic signal detection (Pagel's lambda, Bloomberg's K)
- Trait profiling across phylogeny
- Interactive phylogenetic tree visualizations
- Complete test suite with pytest
- Docker container support
- Google Colab notebook
- Comprehensive documentation
- CI/CD workflows
- Pre-commit hooks
- Example datasets with Newick trees

### Changed
- Optimized GitHub Actions workflows to reduce runner minutes
- Docker builds now only run on releases and manual triggers
- Updated mypy configuration to Python 3.9

### Fixed
- Fixed 12 code quality issues identified by ruff linting
- Fixed class naming inconsistency (PhyloAnalyzer → PhyloTraitAnalyzer)
- Fixed ambiguous variable name 'l' → 'label'
- Fixed all bare except clauses with specific exception handling
- Fixed type annotation issues for mypy compliance
- Fixed import statements in __init__.py and cli.py
- Fixed example data and tree file inclusion

### Features
- Phylogenetic tree parsing and manipulation
- Ancestral state reconstruction
- Phylogenetic conservation analysis
- Tree-based clustering methods
- Statistical testing with phylogenetic correction
- Interactive HTML reports with phylogenetic visualizations

### Technical Details
- Python 3.8+ support
- BioPython for phylogenetic tree handling
- Reproducible analyses with fixed random seeds
- Command-line interface and Python API
- Docker containerization

## Project History

This tool was developed as part of the StrepSuis Suite for bacterial genomics research,
with a focus on *Streptococcus suis* but applicable to any bacterial species.

[Unreleased]: https://github.com/MK-vet/strepsuis-phylotrait/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/MK-vet/strepsuis-phylotrait/releases/tag/v1.0.0
