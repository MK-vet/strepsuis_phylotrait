# Performance Validation Report - strepsuis-phylotrait

**Generated:** 2026-01-29T07:11:54.825302

---

## Benchmark Summary

| Operation | Data Size | Time (ms) | Throughput |
|-----------|-----------|-----------|------------|
| Jaccard Distance (n=50) | 1225 pairs | 0.1 | 11835745 pairs/s |
| Jaccard Distance (n=100) | 4950 pairs | 0.17 | 28530262 pairs/s |
| Jaccard Distance (n=200) | 19900 pairs | 0.74 | 27034370 pairs/s |
| Jaccard Distance (n=500) | 124750 pairs | 3.79 | 32939903 pairs/s |
| Hamming Distance (n=50) | 1225 pairs | 0.05 | 27041945 pairs/s |
| Hamming Distance (n=100) | 4950 pairs | 0.05 | 95192329 pairs/s |
| Hamming Distance (n=200) | 19900 pairs | 0.34 | 59279119 pairs/s |
| Hamming Distance (n=500) | 124750 pairs | 0.78 | 159323123 pairs/s |
| Hierarchical Clustering (n=50) | 50 samples | 0.23 | 215983 samples/s |
| Hierarchical Clustering (n=100) | 100 samples | 0.2 | 508906 samples/s |
| Hierarchical Clustering (n=200) | 200 samples | 0.33 | 611247 samples/s |
| Hierarchical Clustering (n=500) | 500 samples | 1.65 | 302188 samples/s |
| Correlation Matrix (f=20) | 20x20 | 0.13 | 3179650 elem/s |
| Correlation Matrix (f=50) | 50x50 | 0.29 | 8704734 elem/s |
| Correlation Matrix (f=100) | 100x100 | 1.2 | 8310479 elem/s |
| Correlation Matrix (f=200) | 200x200 | 4.29 | 9332711 elem/s |
| Richness Calculation (n=100) | 100 samples | 0.06 | 1785714 samples/s |
| Richness Calculation (n=500) | 500 samples | 0.04 | 12165461 samples/s |
| Richness Calculation (n=1000) | 1000 samples | 0.15 | 6811989 samples/s |
| Richness Calculation (n=2000) | 2000 samples | 0.05 | 38167903 samples/s |
| Beta Diversity (n=50) | 1225 pairs | 4.65 | 263662 pairs/s |
| Beta Diversity (n=100) | 4950 pairs | 18.21 | 271778 pairs/s |
| Beta Diversity (n=200) | 19900 pairs | 73.18 | 271944 pairs/s |
| Permutation Test (n=100) | 100 permutations | 2.06 | 48428 perm/s |
| Permutation Test (n=500) | 500 permutations | 8.65 | 57772 perm/s |
| Permutation Test (n=1000) | 1000 permutations | 18.28 | 54699 perm/s |

---

## Scalability Analysis

### Jaccard Distance

| Data Size | Time (ms) |
|-----------|----------|
| 50 | 0.1 |
| 100 | 0.17 |
| 200 | 0.74 |
| 500 | 3.79 |

### Hamming Distance

| Data Size | Time (ms) |
|-----------|----------|
| 50 | 0.05 |
| 100 | 0.05 |
| 200 | 0.34 |
| 500 | 0.78 |

### Hierarchical Clustering

| Data Size | Time (ms) |
|-----------|----------|
| 50 | 0.23 |
| 100 | 0.2 |
| 200 | 0.33 |
| 500 | 1.65 |

### Correlation Matrix

| Data Size | Time (ms) |
|-----------|----------|
| 20 | 0.13 |
| 50 | 0.29 |
| 100 | 1.2 |
| 200 | 4.29 |

### Richness Calculation

| Data Size | Time (ms) |
|-----------|----------|
| 100 | 0.06 |
| 500 | 0.04 |
| 1000 | 0.15 |
| 2000 | 0.05 |

### Beta Diversity

| Data Size | Time (ms) |
|-----------|----------|
| 50 | 4.65 |
| 100 | 18.21 |
| 200 | 73.18 |

### Permutation Test

| Data Size | Time (ms) |
|-----------|----------|
| 100 | 2.06 |
| 500 | 8.65 |
| 1000 | 18.28 |

