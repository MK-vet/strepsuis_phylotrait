# Synthetic Data Generation Methodology for Phylogenetic Analysis

## Overview

This document describes the statistical methodology used to generate synthetic
data for validating phylogenetic clustering and trait analysis.

## Generation Parameters

- **Number of strains**: 200
- **Number of clusters**: 4
- **MIC features**: 10
- **AMR gene features**: 15
- **Virulence factor features**: 10
- **Trait heritability**: 0.60
- **Tree depth**: 0.50
- **Random seed**: 42

## Statistical Methods Used

### 1. Tree Generation

A random phylogenetic tree is generated with:
- Clustered topology: strains within the same cluster share a more recent common ancestor
- Exponential branch lengths: shorter branches within clusters, longer between clusters
- Birth-death process: tree structure follows coalescent-like dynamics

### 2. Trait Evolution (Heritability Model)

Binary traits are generated with phylogenetic signal:
- Each cluster has characteristic prevalence patterns
- Heritability parameter (0.60) controls how strongly traits follow phylogeny
- Higher heritability = more phylogenetic signal

Formula:
  P(trait=1|cluster) = heritability × cluster_pattern + (1-heritability) × base_prevalence

### 3. Cluster Structure

True clusters are embedded in the tree structure:
- Each cluster forms a monophyletic clade
- Cluster boundaries correspond to deep splits in the tree

## Ground Truth

### Cluster Distribution
- Cluster 1: 50 strains (25.0%)
- Cluster 2: 50 strains (25.0%)
- Cluster 3: 50 strains (25.0%)
- Cluster 4: 50 strains (25.0%)

## Expected Analysis Performance

- Tree-aware clustering should recover ~60% of true clusters
- Phylogenetic diversity should reflect tree structure
- Trait associations should show phylogenetic signal

## References

1. Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
2. Pagel, M. (1999). Inferring the historical patterns of biological evolution.

## Generation Timestamp

2026-01-20T00:27:38.081372

---
*This data was generated for validation and testing purposes only.*
