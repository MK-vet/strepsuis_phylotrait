# Performance Benchmarks Report - strepsuis-phylotrait

**Generated:** 2026-01-29T07:11:54.826262
**Total Benchmarks:** 25

---

## Benchmark Results

| Operation | Samples | Features | Time (s) | Throughput (samples/s) |
|-----------|---------|----------|----------|------------------------|
| Tree Parsing | 10 | 1 | 0.000 | 167103.7 |
| Tree Parsing | 50 | 1 | 0.000 | 356053.0 |
| Tree Parsing | 100 | 1 | 0.000 | 411206.3 |
| Distance Matrix | 10 | 10 | 0.001 | 12125.8 |
| Distance Matrix | 30 | 30 | 0.016 | 1855.6 |
| Distance Matrix | 50 | 50 | 0.069 | 723.9 |
| Faith's PD | 10 | 1 | 0.000 | 257319.3 |
| Faith's PD | 30 | 1 | 0.000 | 382459.3 |
| Faith's PD | 50 | 1 | 0.000 | 394943.9 |
| MPD/MNTD | 10 | 10 | 0.001 | 10546.4 |
| MPD/MNTD | 30 | 30 | 0.016 | 1882.1 |
| MPD/MNTD | 50 | 50 | 0.233 | 214.5 |
| Hierarchical Clustering | 10 | 10 | 0.001 | 7623.2 |
| Hierarchical Clustering | 30 | 30 | 0.017 | 1742.7 |
| Hierarchical Clustering | 50 | 50 | 0.111 | 449.2 |
| Permutation Test | 50 | 100 | 0.002 | 24684.0 |
| Permutation Test | 100 | 200 | 0.004 | 27209.2 |
| Bootstrap CI | 50 | 500 | 0.004 | 12391.6 |
| Bootstrap CI | 100 | 500 | 0.004 | 22959.8 |
| Bootstrap CI | 200 | 500 | 0.004 | 44976.7 |
| Correlation Matrix | 20 | 10 | 0.000 | 198312.2 |
| Correlation Matrix | 50 | 20 | 0.000 | 788403.0 |
| Correlation Matrix | 100 | 30 | 0.000 | 1503334.8 |
| Full Pipeline | 20 | 10 | 0.011 | 1761.9 |
| Full Pipeline | 50 | 10 | 0.145 | 345.8 |

---

## Performance Summary

### Tree Parsing

- **Average Throughput:** 311454.3 samples/s
- **Scalability:** Tested with 10-100 samples

### Distance Matrix

- **Average Throughput:** 4901.8 samples/s
- **Scalability:** Tested with 10-50 samples

### Faith's PD

- **Average Throughput:** 344907.5 samples/s
- **Scalability:** Tested with 10-50 samples

### MPD/MNTD

- **Average Throughput:** 4214.3 samples/s
- **Scalability:** Tested with 10-50 samples

### Hierarchical Clustering

- **Average Throughput:** 3271.7 samples/s
- **Scalability:** Tested with 10-50 samples

### Permutation Test

- **Average Throughput:** 25946.6 samples/s
- **Scalability:** Tested with 50-100 samples

### Bootstrap CI

- **Average Throughput:** 26776.0 samples/s
- **Scalability:** Tested with 50-200 samples

### Correlation Matrix

- **Average Throughput:** 830016.7 samples/s
- **Scalability:** Tested with 20-100 samples

### Full Pipeline

- **Average Throughput:** 1053.9 samples/s
- **Scalability:** Tested with 20-50 samples

