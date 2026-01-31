# Real Data Validation Report - strepsuis-phylotrait

**Generated:** 2026-01-29T07:11:54.824505
**Data Source:** S. suis strains (Snp_tree.newick, AMR_genes.csv, Virulence.csv)
**Total Tests:** 9
**Passed:** 9
**Coverage:** 100.0%

---

## Statistical Validation Results

| Test | Expected | Actual | Status |
|------|----------|--------|--------|
| Tree Loading | >0 terminals | 91 terminals | ✅ PASS |
| Tree Structure | Valid structure | 91 tips, 81 internal | ✅ PASS |
| Distance Matrix Properties | Symmetric, non-negative, zero diagonal | Sym=True, Non-neg=True, Diag0=True | ✅ PASS |
| Trait Prevalence | [0%, 100%] | Range: [1.1%, 60.4%] | ✅ PASS |
| Correlation Matrix Properties | Symmetric, diagonal=1, [-1,1] | Sym=True, Diag1=True, Range=True | ✅ PASS |
| Trait-Tree Matching | ≥10 common strains | 91 common | ✅ PASS |
| Trait Diversity | Valid statistics | Mean richness: 2.5±2.0 | ✅ PASS |
| Beta Diversity | [0, 1] | Mean Jaccard: 0.664 | ✅ PASS |
| Trait Variability | ≥1 variable trait | 5/5 variable | ✅ PASS |

---

## Biological Validation Results

### Phylogenetic Tree

**Description:** SNP-based phylogenetic tree of S. suis strains

**Result:** 91 terminal nodes (strains)

**Interpretation:** Tree represents evolutionary relationships based on core genome SNPs.

### Phylogenetic Distances

**Description:** Pairwise phylogenetic distances between strains

**Result:** Mean distance: 0.1270

**Interpretation:** Larger distances indicate more divergent strains. Clustering of distances may indicate population structure.

### AMR Gene Prevalence

**Description:** Prevalence of AMR genes across strains

**Result:** Top 3: erm(B):60%, tet(O):44%, ant(6)-Ia:30%

**Interpretation:** High prevalence genes may be core resistance determinants or located on successful mobile elements.

### Gene Correlation

**Description:** Strongest correlation between AMR genes

**Result:** lnu(B) vs lsa(E): r=1.000

**Interpretation:** Strong positive correlation suggests co-selection or genetic linkage. Negative correlation may indicate functional redundancy.

### Phylogenetic Signal

**Description:** Matching between phylogenetic tree and trait data

**Result:** 91 strains with both tree and trait data

**Interpretation:** Phylogenetic signal analysis requires matching strain identifiers between tree and trait data.

### AMR Gene Richness

**Description:** Number of AMR genes per strain

**Result:** Mean: 2.5 ± 2.0, Range: 0-9

**Interpretation:** Higher richness indicates more extensive resistance profiles. Variation suggests diverse selection pressures.

### Beta Diversity

**Description:** Average Jaccard distance between strain AMR profiles

**Result:** Mean Jaccard distance: 0.664

**Interpretation:** Higher values indicate more diverse AMR profiles. Lower values suggest homogeneous population.

### Trait Evolution Potential

**Description:** Traits with sufficient variability for evolutionary analysis

**Result:** 5 traits show variation across strains

**Interpretation:** Variable traits can be analyzed for phylogenetic signal and evolution rate.

