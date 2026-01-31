# StrepSuis-PhyloTrait Tutorial

## Quick Start Guide

This tutorial will guide you through using StrepSuis-PhyloTrait for phylogenetic trait analysis.

## Table of Contents

1. [Installation](#installation)
2. [Data Preparation](#data-preparation)
3. [Running Analysis](#running-analysis)
4. [Understanding Results](#understanding-results)
5. [Advanced Usage](#advanced-usage)

## Installation

### Option 1: pip install (recommended)

```bash
pip install strepsuis-phylotrait
```

### Option 2: From source

```bash
git clone https://github.com/MK-vet/MKrep.git
cd MKrep/separated_repos/strepsuis-phylotrait
pip install -e .
```

### Verify installation

```python
import strepsuis_phylotrait
print(strepsuis_phylotrait.__version__)
```

## Data Preparation

### Required Files

Your data directory should contain:

1. **Snp_tree.newick** - Phylogenetic tree in Newick format
2. **AMR_genes.csv** - Antimicrobial resistance genes
3. **Virulence.csv** - Virulence factors

### File Formats

#### Newick Tree
```
((A:1.0,B:1.0):0.5,(C:1.0,D:1.0):0.5);
```

#### CSV Files
- First column: `Strain_ID` (must match tree tip labels)
- Binary values: 0 (absent) or 1 (present)
- UTF-8 encoding

Example:
```csv
Strain_ID,Gene1,Gene2,Gene3
A,1,0,1
B,0,1,1
C,1,1,0
D,0,0,1
```

## Running Analysis

### Command Line Interface

```bash
# Basic usage
strepsuis-phylotrait --data-dir ./data --output ./results

# With custom parameters
strepsuis-phylotrait \
  --data-dir ./data \
  --output ./results \
  --n-clusters 4 \
  --bootstrap 1000
```

### Python API

```python
from strepsuis_phylotrait import PhyloTraitAnalyzer

# Initialize
analyzer = PhyloTraitAnalyzer(
    data_dir="./data",
    output_dir="./results"
)

# Run analysis
results = analyzer.run_complete_analysis()

# Generate reports
analyzer.generate_html_report(results)
analyzer.generate_excel_report(results)
```

## Understanding Results

### Output Files

1. **HTML Report** - Interactive phylogenetic visualizations
2. **Excel Report** - Trait statistics and cluster assignments
3. **PNG Charts** - Tree plots, UMAP projections

### Key Metrics

- **Faith's PD**: Phylogenetic diversity
- **MPD**: Mean pairwise distance
- **MNTD**: Mean nearest taxon distance
- **Fritz-Purvis D**: Phylogenetic signal

## Advanced Usage

### Custom Configuration

```python
from strepsuis_phylotrait import Config, PhyloTraitAnalyzer

config = Config(
    n_clusters=4,
    bootstrap_iterations=1000,
    fdr_alpha=0.05
)

analyzer = PhyloTraitAnalyzer(
    data_dir="./data",
    output_dir="./results",
    config=config
)
```

### Using Innovations

#### Trait Evolution Rate Estimation

```python
from strepsuis_phylotrait.phylo_analysis_core import estimate_evolution_rate

# Estimate rate for a trait
rate = estimate_evolution_rate(tree, trait_values)
print(f"Evolution rate: {rate}")
```

#### Phylogenetic Trait Correlation

```python
from strepsuis_phylotrait.phylo_analysis_core import phylogenetic_correlation

# Calculate phylogenetically-corrected correlation
corr, p_value = phylogenetic_correlation(
    tree, trait1, trait2
)
```

### Tree-Aware Clustering

```python
from Bio import Phylo
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# Load tree
tree = Phylo.read("tree.newick", "newick")

# Get distance matrix
taxa = [t.name for t in tree.get_terminals()]
dist_matrix = [[tree.distance(t1, t2) for t2 in taxa] for t1 in taxa]

# Cluster
Z = linkage(squareform(dist_matrix), method='average')
clusters = fcluster(Z, t=4, criterion='maxclust')
```

## Troubleshooting

### Common Issues

1. **Tree/Data Mismatch**: Ensure Strain_IDs match tree tip labels
2. **Invalid Newick**: Check tree format and branch lengths
3. **Memory Error**: Use subset of taxa for large trees

### Getting Help

- GitHub Issues: https://github.com/MK-vet/strepsuis-phylotrait/issues
- Documentation: See README.md and USER_GUIDE.md

## Performance Tips

Based on our benchmarks:

| Operation | Throughput |
|-----------|------------|
| Tree Parsing | 147,262 samples/s |
| Distance Matrix | 2,338 samples/s |
| Faith's PD | 170,009 samples/s |
| Full Pipeline | 1,027 samples/s |

## Phylogenetic Signal Interpretation

### Fritz-Purvis D Statistic

- **D = 0**: Brownian motion (strong phylogenetic signal)
- **D = 1**: Random distribution (no phylogenetic signal)
- **D < 0**: More clustered than Brownian
- **D > 1**: Overdispersed

### Practical Interpretation

| D Value | Interpretation |
|---------|----------------|
| < 0.3 | Strong phylogenetic signal |
| 0.3 - 0.7 | Moderate signal |
| > 0.7 | Weak/no signal |

## Next Steps

- Read [ALGORITHMS.md](ALGORITHMS.md) for details on novel features and algorithms
- See [VALIDATION.md](VALIDATION.md) for statistical validation
- Check [BENCHMARKS.md](BENCHMARKS.md) for performance data
