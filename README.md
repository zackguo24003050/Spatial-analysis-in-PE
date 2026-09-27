# Spatial Analysis at the Maternal-Fetal Interface in Severe Preeclampsia

This repository contains the analysis workflow for a pilot Xenium spatial
transcriptomics project investigating trophoblast-endothelial organization and
cell-cell communication in severe preeclampsia (PE).

The project began with an initial comparison of one control and one disease
sample, and then expanded to a six-donor dataset containing two controls and
four PE samples. The final single-donor and multisample SpatialCellChat
calculations were run as R scripts
on the Myriad computing cluster. The R Markdown files contain
detailed script and interpretation. Major results are shown below in 
"Analysis Workflow and Preliminary Results".

**Supervisor:** Yara Elana Sanchez Corrales<br>
**Principal investigator:** Sergi Castellano Hereza<br>
**Analysis status:** preliminary study; additional samples and refined cell-type
annotations are in progress<br>
**Main question:** are extravillous trophoblast-endothelial spatial relationships
and candidate ligand-receptor communication programs altered in severe PE?

## Project Overview

Preeclampsia is associated with abnormal placentation and impaired spiral artery
remodeling. During normal pregnancy, extravillous trophoblasts (EVTs) invade the
maternal decidua and myometrium and interact with vascular endothelial cells as
spiral arteries are converted into low-resistance vessels. Spatial
transcriptomics makes it possible to examine both the physical organization of
these cells and the signaling programs that may operate between them.

The analysis used a targeted 479-gene Xenium panel. The initial analysis
quantified nearest neighbors and intercellular distances while comparing one
control donor (FVB) with one disease donor (FAM). Ligand-receptor analysis was
explored with LIANA, including a nearest-neighbor-restricted version. Because
that restriction did not produce sufficiently useful spatial results, the final
workflow used SpatialCellChat to incorporate cell coordinates directly into
communication inference.

The cohort was subsequently expanded to six donors:

| Condition | Donors | Number of donors |
|---|---|---:|
| Control | FVB, FVQ | 2 |
| Severe PE | FAM, FCM, FEP, FVS | 4 |

The current results should be treated as hypothesis-generating. They establish
a working analysis framework and prioritize candidate pathways for a larger
cohort planned for future study.

```text
01 Xenium spatial data and cell-type annotations
                         |
                         v
02 Nearest-neighbor preference and distance analysis
                         |
                         v
03 Exploratory LIANA ligand-receptor analysis
   whole tissue and k-nearest-neighbor-restricted subsets
                         |
                         v
04 Donor and condition level comparison,
   and candidate ligand-receptor checks
```

## Data and Analysis Design

Contact-dependent signaling uses a 20 micrometer contact range. Secreted
signaling uses a 250 micrometer interaction range. Both workflows calculate
cell-group pathway probabilities and individual-cell pathway probabilities.
The CellChat human database is supplemented with the candidate
`CLEC11A-KIT` interaction described in the maternal-fetal interface literature
[2].

The original data are not intended to be shown before formal result is published.

## Analysis Workflow and Preliminary Results

### 01. Spatial neighborhood exploration

[`Neigbour_analysis_exploration.Rmd`](Neigbour_analysis_exploration.Rmd)
compares endothelial extravillous trophoblasts (eEVTs) and maternal endothelial
cells in one control donor (FVB) and one PE donor (FAM).

**Nearest-neighbor preference (NEP).** Within each donor and tissue region
(decidua or muscle), I used cell-centroid coordinates to find the *n* nearest
other cells for every eEVT and endothelial cell. I counted directed pairs from
each query cell type to each neighbor cell type, then divided by all neighbor
pairs from that query type in the same region. For example, the endothelial to
eEVT percentage is the number of endothelial cells whose nearest neighbor is
an eEVT divided by the number of endothelial query cells when *n* = 1. The
heatmaps below use *n* = 1; each row sums to 100%. 

![Nearest-neighbor cell-type proportions for FVB and FAM, shown separately for decidua and muscle](figures/nearest-neighbor-proportions.png)

*Rows show the query cell type and columns show the neighbor cell type. In the
decidua, 6.7% of endothelial cells in FVB had an eEVT as their nearest neighbor,
compared with 0.4% in FAM. These proportions describe the observed neighborhood
composition and depend on how many eEVTs are present.*

**COZI.** I also ran conditional neighborhood preference analysis separately
for FVB and FAM using each cell's three nearest neighbors and 300 label
permutations. This calculation uses each donor as a whole rather than separating
decidua and muscle. In the dot plots, dot size is the conditional cell ratio (the
fraction of source cells with at least one neighbor of the indicated type).
Color is the z-score for the conditional neighbor count relative to the
permuted cell labels: positive values indicate a higher value than expected
under that null model, and negative values indicate a lower value.

![COZI conditional cell ratios and permutation z-scores for FVB and FAM](figures/cozi-scores.png)

*The initial pair suggests a weaker permutation-standardized eEVT self-neighbor
signal in FAM.
Because this comparison includes only one donor per condition, the figure does
not establish a reproducible PE effect. The marked difference in eEVT abundance
also complicates interpretation of endothelial to eEVT proportions.*

### 02. Exploratory LIANA analysis

[`Ligand-Receptor Exploration organized.Rmd`](Ligand-Receptor%20Exploration%20organized.Rmd)
contains the LIANA workflow, including eEVT-endothelial-focused and
nearest-neighbor-restricted analyses.

**Result:** LIANA prioritized NOTCH- and LEP-related candidates in
the initial control-disease comparison. Restricting the expression data to very
small spatial neighborhoods did not yield sufficiently stable or interpretable
results, so this branch was retained as method exploration rather than the
final communication analysis.

### 03. Initial SpatialCellChat analysis

The initial FVB-FAM analysis established the SpatialCellChat workflow. It
separated cell-cell contact interactions from secreted signaling and calculated
ligand-receptor and pathway probabilities using Xenium expression and cell
coordinates.

**Preliminary result:** the initial pair suggested relatively stronger NOTCH
and CDH5-associated communication in the control sample and relatively stronger
inflammatory cytokine signaling in the disease sample. This analysis served as
method development and motivated expansion to all six donors.

### 04. Six-donor single-sample calculation

The per-donor calculation was performed on the Myriad cluster with
[`SpatialCellChat_single_sample_compute_myriad.R`](SpatialCellChat_single_sample_compute_myriad.R).
The resulting independent donor objects supported donor-level pathway
summaries, communication heatmaps, hotspot visualization, and individual-cell
pathway maps.

### 05. Multisample condition calculation

The condition-level calculation was performed on Myriad with
[`SpatialCellChat_multisample_compute_myriad.R`](SpatialCellChat_multisample_compute_myriad.R).
The script groups FVB and FVQ as control and FAM, FCM, FEP, and FVS as disease.
[`Comparison-by-Cellchat-multiple-sample.Rmd`](Comparison-by-Cellchat-multiple-sample.Rmd)
contains the downstream control-disease comparison, including donor-averaged
interaction-count summaries and pooled condition-level comparisons.

**Result:** the six-donor analysis prioritized NOTCH, NCAM, and VEGF
as control-associated candidate programs, while the disease condition showed
stronger inflammatory signaling candidates, including TNF-, IL1-, IL2-, and
IL4-related pathways. Individual-cell pathway maps were used to examine whether
these scores localized to trophoblast or endothelial compartments.

## Interpretation of the Multisample Comparison

The repository contains two related but distinct forms of comparison:

- **Donor-level summaries** are calculated from the six independently inferred
  SpatialCellChat objects.
- **Pooled condition objects** retain donor identity through the `samples`
  field but infer one communication network for all control cells and one for
  all disease cells. RankNet, centrality, and several differential network plots
  describe these pooled objects.

The pooled plots are useful for pathway prioritization, but they do not
constitute a donor-level statistical test. The unequal cohort size of two
controls and four disease donors should also be considered when interpreting
condition-level differences.

## Overall Conclusion

This pilot project established an end-to-end workflow for spatial neighborhood
analysis and ligand-receptor inference in Xenium data from the maternal-fetal
interface. The work progressed from an initial one-control, one-disease
comparison to a six-donor SpatialCellChat analysis.

Across the exploratory analyses, NOTCH, NCAM, CLEA11A, VEGF, and inflammatory signaling
emerged as candidates for follow-up. The evidence remains preliminary, but it
provides a focused set of hypotheses and a reusable computational framework for
the group's planned expanded study.

## Limitations and Next Steps

- The cohort contains only two control and four disease donors.
- A targeted 479-gene panel limits the number of ligand-receptor pairs that can
  be detected.
- Current cell-type annotations remain incomplete and may combine biologically
  distinct trophoblast, endothelial, stromal, or immune subpopulations.
- Additional donors, expanded annotations, donor-level statistical analysis,
  and independent biological validation are required before drawing firm
  disease-mechanism conclusions.

## Repository Structure

```text
.
+-- Neigbour_analysis_exploration.Rmd
+-- Ligand-Receptor Exploration organized.Rmd
+-- Comparison-by-Cellchat-multiple-sample.Rmd
+-- SpatialCellChat_single_sample_compute_myriad.R
+-- SpatialCellChat_multisample_compute_myriad.R
+-- figures/
|   +-- nearest-neighbor-proportions.png
|   +-- cozi-scores.png
+-- README.md
```

## Main Analysis Files

| Stage | File | Purpose |
|---:|---|---|
| 01 | [Neigbour_analysis_exploration.Rmd](Neigbour_analysis_exploration.Rmd) | Nearest-neighbor composition, distance, k sensitivity, and COZI exploration |
| 02 | [Ligand-Receptor Exploration organized.Rmd](Ligand-Receptor%20Exploration%20organized.Rmd) | eEVT-endothelial LIANA workflow and k-nearest-neighbor-restricted exploration |
| 03 | [SpatialCellChat_single_sample_compute_myriad.R](SpatialCellChat_single_sample_compute_myriad.R) | SpatialCellChat calculation for each donor independently |
| 04 | [SpatialCellChat_multisample_compute_myriad.R](SpatialCellChat_multisample_compute_myriad.R) | SpatialCellChat calculation for pooled control and disease objects |
| 05 | [Comparison-by-Cellchat-multiple-sample.Rmd](Comparison-by-Cellchat-multiple-sample.Rmd) | Donor-aware summaries and pooled control-disease network comparisons |

## Software

The analysis is written in R. Packages used across the workflows include:

- Seurat
- SpatialCellChat and CellChat
- LIANA
- RANN and coziR
- Matrix, dplyr, tidyr, and tibble
- ggplot2, patchwork, and ComplexHeatmap
- future

## References

1. Dimitrov D, Türei D, Garrido-Rodriguez M, et al. *Comparison of methods and
   resources for cell-cell communication inference from single-cell RNA-Seq
   data.* Nature Communications. 2022;13:3224.
   [doi:10.1038/s41467-022-30755-0](https://doi.org/10.1038/s41467-022-30755-0)
2. Greenbaum S, Averbukh I, Soon E, et al. *A spatially resolved timeline of the
   human maternal-fetal interface.* Nature. 2023;619:595-605.
   [doi:10.1038/s41586-023-06298-9](https://doi.org/10.1038/s41586-023-06298-9)
3. Jin S, Plikus MV, Nie Q. *CellChat for systematic analysis of cell-cell
   communication from single-cell transcriptomics.* Nature Protocols.
   2025;20:180-219.
   [doi:10.1038/s41596-024-01045-4](https://doi.org/10.1038/s41596-024-01045-4)
4. Schiller C, Ibarra-Arellano MA, Bestak K, Tanevski J, Schapiro D.
   *Comparison and optimization of cellular neighbor preference methods for
   quantitative tissue analysis.* Nature Communications. 2026;17.
   [doi:10.1038/s41467-026-71699-z](https://doi.org/10.1038/s41467-026-71699-z)
