# Putative Disease Gene Identification and Drug Repurposing for Depressive Neurosis

This repository contains the code and documentation for our Master's project in Bioinformatics and Network Medicine (Data Science, a.y. 2024/2025). The project leverages network medicine to identify novel disease genes for depressive neurosis and explores drug repurposing opportunities based on protein-protein interaction (PPI) data and gene-disease associations (GDAs).

## Table of Contents

- [Overview](#overview)
- [Abstract](#abstract)
- [Introduction](#introduction)
- [Materials and Methods](#materials-and-methods)
  - [PPI Data and Gene-Disease Associations](#ppi-data-and-gene-disease-associations)
  - [Disease Subnetwork Construction](#disease-subnetwork-construction)
  - [Algorithm Comparison and Cross-Validation](#algorithm-comparison-and-cross-validation)
  - [Enrichment Analysis](#enrichment-analysis)
  - [Drug Repurposing](#drug-repurposing)
- [Results and Discussion](#results-and-discussion)
- [Author Contributions](#author-contributions)
- [Programming Languages and Tools](#programming-languages-and-tools)
- [References](#references)
- [License](#license)

## Overview

This project aims to identify putative disease genes associated with depressive neurosis and to repurpose existing drugs for potential therapeutic applications. By constructing a human interactome from BioGRID and integrating curated GDAs, we isolate a disease-specific subnetwork. We then compare multiple gene identification algorithms—DIAMOnD, DiaBLE, and a diffusion-based method—to prioritize candidate genes, followed by enrichment and drug repurposing analyses.

## Abstract

Depressive neurosis is a chronic mood disorder that severely impacts quality of life. Using a network-based approach, we mapped 271 depressive neurosis-related genes onto a human PPI network (filtered to 263 mapped genes). The largest connected component (LCC) of the resulting disease network comprises 197 nodes. We evaluated three algorithms for disease gene identification and found that a diffusion-based approach (with a diffusion time of t = 0.002) outperformed others based on 5-fold cross-validation metrics. Subsequent enrichment analysis and drug repurposing identified TAMOXIFEN and CISPLATIN as promising therapeutic candidates.

## Introduction

Depressive neurosis is characterized by persistent feelings of sadness and hopelessness. Network medicine provides a framework to explore the molecular basis of such complex disorders. In this project, we:
- Construct a comprehensive human interactome using BioGRID data.
- Map curated depressive neurosis GDAs onto the interactome.
- Analyze the resulting disease subnetwork to identify key genes.
- Compare algorithmic approaches for disease gene prioritization.
- Perform enrichment analysis and explore drug repurposing opportunities.

## Materials and Methods

### PPI Data and Gene-Disease Associations

- **PPI Data:**  
  Obtained from the latest BioGRID release, filtering for physical interactions in Homo sapiens (taxon ID 9606) and removing redundant edges and self-loops.

- **Gene-Disease Associations (GDAs):**  
  Curated GDAs for depressive neurosis yielded 271 genes. After verification against the HGNC database, 263 genes were successfully mapped to the interactome.

### Disease Subnetwork Construction

- The 263 mapped genes form the basis of a disease-specific subnetwork.
- The Largest Connected Component (LCC) of this network contains 197 nodes.
- Topological analyses (degree, betweenness, eigenvector, and closeness centrality) highlighted key genes such as **APP**, **STIP1**, and **AKT1**.

### Algorithm Comparison and Cross-Validation

- **Algorithms Evaluated:**
  - **DIAMOnD** and **DiaBLE:** Hypergeometric-based iterative methods.
  - **Diffusion-based Method:** A random walk with restart approach.
  
- **Cross-Validation Framework:**
  - A 5-fold cross-validation was performed. Seed genes were assigned a heat value of 1, while other genes were initialized to 0.
  - Performance metrics (precision, recall, F1-score, and accuracy) were computed on the top 50 predictions for each fold.
  - The diffusion-based method (with t = 0.002) achieved an average F1-score of 0.00777, outperforming the other approaches.

### Enrichment Analysis

- Two gene lists were analyzed:
  - The original set of depressive neurosis genes.
  - The top 100 putative disease genes predicted by the diffusion algorithm.
- Enrichment analysis was performed using Enrichr across multiple categories (GO Biological Process, Molecular Function, Cellular Component, Reactome, and KEGG pathways).

### Drug Repurposing

- The top 20 putative disease genes were cross-referenced with the DGIdb dataset.
- Drugs were ranked based on the number of target genes.
- **Key Candidates:**
  - **TAMOXIFEN** and **CISPLATIN:** Each targeting three genes.
  - **DIETHYLSTILBESTROL:** Also a top candidate but lacking clinical trial support.
- ClinicalTrials.gov searches confirmed TAMOXIFEN (85 trials) and CISPLATIN (21 trials) as established candidates for repurposing, suggesting greater therapeutic potential.

## Results and Discussion

- **Network Analysis:**  
  The LCC of the disease network (197 nodes) provided a robust framework for centrality analyses.
  
- **Algorithm Evaluation:**  
  Despite low absolute performance metrics, the diffusion-based algorithm consistently outperformed DIAMOnD and DiaBLE in retrieving withheld disease genes.
  
- **Enrichment and Drug Repurposing:**  
  Enrichment analysis of the original gene set revealed significant pathways, while drug repurposing identified TAMOXIFEN and CISPLATIN as promising candidates, validated by clinical trial evidence.

## Author Contributions

- **Mayis Atayev (2104359):**  
  Implemented algorithms, conducted computational validations, performed gene identification, enrichment analysis, drug repurposing, and drafted the report.
  
- **Dila Aslan (2113310):**  
  Managed data collection, algorithm implementation, cross-validation, gene identification, drug repurposing, and contributed to report drafting and review.

## Programming Languages and Tools

- **Language:** Python
- **Key Libraries/Tools:**
  - **NetworkX:** For network construction and analysis.
  - **Pandas and NumPy:** For data manipulation.
  - **Scikit-learn:** For performance metric evaluation.
  - **Matplotlib and Seaborn:** For visualizations.
  - **Enrichr:** For functional enrichment analysis.
  - **DGIdb:** For drug-gene interaction analysis.
  - **ClinicalTrials.gov:** For clinical trial validation.

## References

1. Jensen, L. J., et al. *DISEASES: Gene-Disease Associations database.* Available at: [https://diseases.jensenlab.org](https://diseases.jensenlab.org). Accessed: January 2025.
2. Stark, C., et al. *BioGRID: A resource for studying protein and genetic interactions.* Available at: [https://thebiogrid.org](https://thebiogrid.org). Accessed: January 2025.
3. Wagner, A. H., et al. *DGIdb: The Drug-Gene Interaction Database.* Available at: [https://www.dgidb.org](https://www.dgidb.org). Accessed: January 2025.
4. Kuleshov, M. V., et al. *Enrichr: A comprehensive gene set enrichment analysis web server.* Available at: [https://maayanlab.cloud/Enrichr](https://maayanlab.cloud/Enrichr). Accessed: January 2025.
5. The Gene Ontology Consortium. *Gene Ontology: Tool for the unification of biology.* Available at: [http://geneontology.org](http://geneontology.org). Accessed: January 2025.
6. Jassal, B., et al. *Reactome: A curated knowledgebase of biological pathways.* Available at: [https://reactome.org](https://reactome.org). Accessed: January 2025.
7. Kanehisa, M., and Goto, S. *KEGG: Kyoto Encyclopedia of Genes and Genomes.* Available at: [https://www.genome.jp/kegg](https://www.genome.jp/kegg). Accessed: January 2025.
8. Clinical Trial on TAMOXIFEN: [https://clinicaltrials.gov/study/NCT00667121?cond=Depressive%20neurosis&term=TAMOXIFEN&rank=1](https://clinicaltrials.gov/study/NCT00667121?cond=Depressive%20neurosis&term=TAMOXIFEN&rank=1). Accessed: January 2025.
9. Clinical Trial on CISPLATIN: [https://clinicaltrials.gov/study/NCT00005850?cond=Depressive%20neurosis&term=CISPLATIN&rank=1](https://clinicaltrials.gov/study/NCT00005850?cond=Depressive%20neurosis&term=CISPLATIN&rank=1). Accessed: January 2025.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
