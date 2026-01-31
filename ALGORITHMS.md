# Algorithms Documentation

This document provides detailed algorithmic descriptions and Big-O complexity analysis
for the key computational methods in the StrepSuis-PhyloTrait module.

## Overview

All algorithms in the StrepSuis-PhyloTrait module are designed with:
- **Reproducibility**: Fixed random seeds for deterministic results
- **Numerical stability**: Careful handling of edge cases and numerical precision
- **Scalability**: Efficient implementations suitable for datasets up to 10,000+ strains

## 1. Faith's Phylogenetic Diversity (PD)

### Purpose
Quantify the phylogenetic diversity of a set of taxa based on tree branch lengths.

### Algorithm
```
FUNCTION faiths_pd(tree, taxa_set):
    INPUT:
        tree: Phylogenetic tree (Newick format)
        taxa_set: Set of taxa to calculate PD for
    
    OUTPUT:
        pd_value: Sum of branch lengths connecting taxa
    
    PROCEDURE:
        # Find minimum spanning tree connecting taxa
        subtree = empty tree
        
        FOR EACH taxon in taxa_set:
            path = path_to_root(tree, taxon)
            FOR EACH edge in path:
                IF edge not in subtree:
                    subtree.add(edge)
        
        # Sum branch lengths
        pd_value = sum(edge.length for edge in subtree)
        
        RETURN pd_value

FUNCTION faiths_pd_vectorized(tree, trait_matrix):
    # Calculate PD for each trait
    pd_values = []
    
    FOR EACH trait in trait_matrix.columns:
        positive_taxa = taxa where trait == 1
        pd = faiths_pd(tree, positive_taxa)
        pd_values.append(pd)
    
    RETURN pd_values
```

### Complexity
- Time: O(n × d) per trait, where n = taxa, d = tree depth
- Space: O(n) for storing visited edges

## 2. Phylogenetic Distance Matrix

### Purpose
Compute pairwise phylogenetic distances between all taxa.

### Algorithm
```
FUNCTION compute_phylo_distance_matrix(tree):
    INPUT:
        tree: Phylogenetic tree
    
    OUTPUT:
        dist_matrix: n × n distance matrix
    
    PROCEDURE:
        n = number of taxa
        dist_matrix = zeros(n, n)
        
        # Pre-compute distances to root for efficiency
        dist_to_root = {}
        lca_depths = {}  # Lowest common ancestor depths
        
        FOR EACH taxon in tree.leaves:
            dist_to_root[taxon] = path_length(taxon, root)
        
        # Compute pairwise distances using LCA
        FOR i = 0 TO n-1:
            FOR j = i+1 TO n-1:
                lca = lowest_common_ancestor(taxa[i], taxa[j])
                d_i_lca = dist_to_root[taxa[i]] - dist_to_root[lca]
                d_j_lca = dist_to_root[taxa[j]] - dist_to_root[lca]
                
                dist_matrix[i, j] = d_i_lca + d_j_lca
                dist_matrix[j, i] = dist_matrix[i, j]
        
        RETURN dist_matrix
```

### Complexity
- Time: O(n² × d) for n taxa and tree depth d
- Space: O(n²) for distance matrix

## 3. Tree-Aware Clustering

### Purpose
Cluster taxa respecting phylogenetic relationships.

### Algorithm
```
FUNCTION phylo_clustering(tree, trait_matrix, k):
    INPUT:
        tree: Phylogenetic tree
        trait_matrix: Binary trait matrix
        k: Number of clusters
    
    OUTPUT:
        cluster_labels: Array of cluster assignments
    
    PROCEDURE:
        # Compute phylogenetic distance matrix
        phylo_dist = compute_phylo_distance_matrix(tree)
        
        # Compute trait distance (Jaccard)
        trait_dist = pdist(trait_matrix, metric='jaccard')
        trait_dist = squareform(trait_dist)
        
        # Combined distance (weighted)
        alpha = 0.5  # Weight for phylogenetic distance
        combined_dist = alpha * phylo_dist + (1 - alpha) * trait_dist
        
        # Hierarchical clustering
        linkage_matrix = linkage(squareform(combined_dist), method='average')
        cluster_labels = cut_tree(linkage_matrix, n_clusters=k)
        
        RETURN cluster_labels.flatten()
```

### Complexity
- Time: O(n² × d + n² log n) for distance computation and clustering
- Space: O(n²)

## 4. Binary Trait Association Analysis

### Purpose
Test for association between binary traits and phylogenetic structure.

### Algorithm
```
FUNCTION phylo_trait_association(tree, trait, n_permutations=999):
    INPUT:
        tree: Phylogenetic tree
        trait: Binary trait vector
        n_permutations: Number of permutations for significance test
    
    OUTPUT:
        observed_pd: Observed phylogenetic diversity
        p_value: Permutation p-value
    
    PROCEDURE:
        # Calculate observed PD for taxa with trait
        positive_taxa = taxa where trait == 1
        observed_pd = faiths_pd(tree, positive_taxa)
        
        # Permutation test
        null_distribution = []
        FOR i = 1 TO n_permutations:
            # Randomly shuffle trait assignments
            shuffled_trait = random_permutation(trait)
            shuffled_positive = taxa where shuffled_trait == 1
            null_pd = faiths_pd(tree, shuffled_positive)
            null_distribution.append(null_pd)
        
        # Calculate p-value (two-tailed)
        n_extreme = count(|null_pd - mean(null)| >= |observed_pd - mean(null)|)
        p_value = (n_extreme + 1) / (n_permutations + 1)
        
        RETURN (observed_pd, p_value)
```

### Complexity
- Time: O(n_permutations × n × d)
- Space: O(n_permutations) for null distribution

## 5. Phylogenetic Signal Metrics

### Purpose
Quantify the degree of phylogenetic signal in trait distributions.

### Algorithms
```
FUNCTION blombergs_k(tree, trait):
    # Measure of phylogenetic signal for continuous traits
    # K = 1 indicates Brownian motion evolution
    # K < 1 indicates less signal than expected
    # K > 1 indicates more signal than expected
    
    # Compute phylogenetic variance-covariance matrix
    C = phylo_vcv(tree)
    
    # Observed mean squared error
    trait_centered = trait - mean(trait)
    mse_observed = trait_centered.T @ inv(C) @ trait_centered / n
    
    # Expected MSE under Brownian motion
    mse_expected = var(trait)
    
    K = (mse_observed / mse_expected) / mean(diag(inv(C)))
    
    RETURN K

FUNCTION consistency_index(tree, trait):
    # CI = min_changes / observed_changes
    # CI = 1 means no homoplasy (trait evolved once)
    
    # Count observed character changes on tree
    observed_changes = parsimony_changes(tree, trait)
    
    # Minimum possible changes (trait states - 1)
    min_changes = len(unique(trait)) - 1
    
    # Maximum possible changes
    max_changes = n - 1
    
    ci = min_changes / observed_changes if observed_changes > 0 else 1.0
    
    RETURN ci

FUNCTION retention_index(tree, trait):
    # RI measures retention of synapomorphy
    
    observed = parsimony_changes(tree, trait)
    min_changes = len(unique(trait)) - 1
    max_changes = n - 1
    
    IF max_changes == min_changes:
        RETURN 1.0
    
    ri = (max_changes - observed) / (max_changes - min_changes)
    
    RETURN ri
```

### Complexity
- Time: O(n² + n × d) for K, O(n × d) for CI/RI
- Space: O(n²) for variance-covariance matrix

## 6. Ancestral State Reconstruction

### Purpose
Infer ancestral trait states at internal nodes.

### Algorithm
```
FUNCTION marginal_ancestral_reconstruction(tree, trait):
    INPUT:
        tree: Phylogenetic tree
        trait: Binary trait vector at tips
    
    OUTPUT:
        ancestral_states: Dictionary of node -> probability of state 1
    
    PROCEDURE:
        # Fitch parsimony for initial states
        # Bottom-up pass
        FOR node in tree.postorder():
            IF node.is_leaf():
                node.state_set = {trait[node.name]}
            ELSE:
                child_states = [child.state_set for child in node.children]
                intersection = set_intersection(child_states)
                IF intersection:
                    node.state_set = intersection
                ELSE:
                    node.state_set = set_union(child_states)
        
        # Top-down pass
        root.final_state = min(root.state_set)
        FOR node in tree.preorder():
            IF not node.is_root():
                parent_state = node.parent.final_state
                IF parent_state in node.state_set:
                    node.final_state = parent_state
                ELSE:
                    node.final_state = min(node.state_set)
        
        # Calculate marginal probabilities using branch lengths
        ancestral_states = {}
        FOR node in tree.nodes():
            ancestral_states[node] = calculate_marginal_prob(node, tree, trait)
        
        RETURN ancestral_states
```

### Complexity
- Time: O(n) for parsimony, O(n²) for marginal reconstruction
- Space: O(n)

## 7. UMAP Embedding for Trait Visualization

### Purpose
Dimensionality reduction for visualizing trait patterns.

### Algorithm
```
FUNCTION umap_embedding(trait_matrix, n_components=2, n_neighbors=15):
    INPUT:
        trait_matrix: Binary matrix of n samples × m traits
        n_components: Output dimensions
        n_neighbors: Number of neighbors for graph
    
    OUTPUT:
        embedding: n × n_components coordinate matrix
    
    PROCEDURE:
        # 1. Construct fuzzy simplicial set (nearest neighbor graph)
        distances = pairwise_distances(trait_matrix, metric='jaccard')
        knn_graph = construct_knn_graph(distances, n_neighbors)
        
        # 2. Compute fuzzy set membership
        sigmas = compute_local_connectivity(knn_graph)
        membership = exp(-(distances - rho) / sigmas)
        
        # 3. Symmetrize the graph
        adjacency = membership + membership.T - membership * membership.T
        
        # 4. Initialize low-dimensional embedding
        embedding = random_init(n_samples, n_components)
        
        # 5. Optimize embedding using SGD
        FOR epoch = 1 TO n_epochs:
            FOR (i, j, weight) in sample_edges(adjacency):
                # Attractive force between connected points
                grad = compute_attractive_gradient(embedding[i], embedding[j], weight)
                embedding[i] += learning_rate * grad
                embedding[j] -= learning_rate * grad
                
                # Repulsive force from random negative sample
                k = random_sample_negative(i, j)
                grad = compute_repulsive_gradient(embedding[i], embedding[k])
                embedding[i] += learning_rate * grad
        
        RETURN embedding
```

### Complexity
- Time: O(n × n_neighbors × n_epochs)
- Space: O(n × n_neighbors) for sparse graph

---

## Scalability Considerations

### Typical Performance

| Operation | 100 strains | 500 strains | 1000 strains |
|-----------|-------------|-------------|--------------|
| Faith's PD (10 traits) | ~0.5s | ~1s | ~2s |
| Distance matrix | ~0.2s | ~2s | ~8s |
| Phylo clustering | ~0.5s | ~5s | ~20s |
| Trait association (100 perms) | ~2s | ~10s | ~30s |
| Full pipeline | ~30s | ~90s | ~180s |

### Memory Optimization

For large datasets (>1000 strains):
1. Use sparse matrix representation for trait data
2. Streaming computation for distance matrices
3. Limit permutation tests to significant traits only

---

## References

1. Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation*, 61(1), 1-10.
2. Blomberg, S. P., et al. (2003). Testing for phylogenetic signal in comparative data. *Evolution*, 57(4), 717-745.
3. Fitch, W. M. (1971). Toward defining the course of evolution. *Systematic Biology*, 20(4), 406-416.
4. McInnes, L., et al. (2018). UMAP: Uniform Manifold Approximation and Projection. *Journal of Open Source Software*, 3(29), 861.
5. Paradis, E., & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. *Bioinformatics*, 35(3), 526-528.

---

**Version**: 1.0.0  
**Last Updated**: 2025-01-15
