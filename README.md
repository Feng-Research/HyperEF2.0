# HyperEF 2.0: Spectral Hypergraph Coarsening via Krylov Subspace Expansion and Resistance-based Local Clustering

**Authors:** Hamed Sajadinia, Zhuo Feng  

---

## Abstract

This paper introduces **HyperEF 2.0**, a scalable framework for **spectral coarsening and clustering** of large-scale hypergraphs through **hyperedge effective resistances**, aiming to decompose hypergraphs into multiple node clusters with a small number of inter-cluster hyperedges. Building on the recent HyperEF framework, our approach offers three primary contributions. Specifically, first, by leveraging the expanded Krylov subspace exploiting both clique and star expansions of hyperedges, we can significantly improve the approximation accuracy of effective resistances. Second, we propose a resistance-based local clustering scheme for merging small isolated nodes into nearby clusters, yielding more balanced clusters with substantially improved conductance. Third, the proposed HyperEF 2.0 enables the integration of resistance-based hyperedge weighting and community detection into a multilevel hypergraph partitioning tool, achieving state-of-the-art performance. Extensive experiments on real-world VLSI benchmarks show that HyperEF 2.0 can more effectively coarsen hypergraphs without compromising their structural properties, while delivering much better solution quality (e.g. conductance) than the state-of-the-art hypergraph coarsening methods, such as HyperEF and HyperSF. Moreover, compared to leading hypergraph partitioners such as hMETIS, SpecPart, MedPart, and KaHyPar, our framework consistently achieves smaller cut sizes. In terms of runtime, HyperEF 2.0 attains up to a $4.5 \times$ speedup over the latest flow-based local clustering algorithm, HyperSF, demonstrating both superior efficiency and partitioning quality.

---

## Requirements

- **Julia Version:** 1.5.3  
- **Required Packages:**
  - `SparseArrays`
  - `LinearAlgebra`
  - `MatrixNetworks v1.0.1`
  - `RandomV06 v0.0.2`

---

## Running HyperEF 2.0

To run HyperEF 2.0 on your own hypergraph:

1. Place your `.hgr` (in hMetis format) hypergraph file in the `data/` directory.
2. In the `src/` directory, Edit `Run_HypperEF_2_0.jl` to include the hypergraph
2. Run:

```bash
julia Run_HypperEF_2_0.jl
```

---

## Output
The output is a clustering file. The i-th line of the file contains the cluster number that the i-th vertex belongs to. Cluster numbers start from 1. 

---

## Citation

If you use this code in your research, please cite the following paper:

Sajadinia, Hamed, and Zhuo Feng. "HyperEF 2.0: Spectral Hypergraph Coarsening via Krylov Subspace Expansion and Resistance-based Local Clustering." Proceedings of the 44st IEEE/ACM International Conference on Computer-Aided Design. 2025.

---

## References

Aghdaei, Ali, and Zhuo Feng. "HyperEF: Spectral hypergraph coarsening by effective-resistance clustering." Proceedings of the 41st IEEE/ACM International Conference on Computer-Aided Design. 2022.

