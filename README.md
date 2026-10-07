---

editor_options: 
  markdown: 
    wrap: 72
---

# R Script for Belowground-drought response analysis

Contact: ***benoit.pichon0\@gmail.com***

<p align="center">
    <img src="https://github.com/bpichon0/Belowground_drought_response/blob/master/Data/idea_project.png" width="800">
</p>

There are 4 numbered scripts. The **0_Functions.R** is some utility funtions.

### Script 1: building the dataset 

- Extracts SPEI (standardized precipitation evapotranspiration index, water deficit) data from NetCDF files for each site (the data is not in the data folders but the link is provided).
- Defines drought/wet thresholds (10th, 25th, 75th, 90th percentiles of the SPEI distribution). This defines extreme, moderate droughts, and moderate and extreme wet events
- Merges SPEI with vegetation data Note that the SPEI is defined as a cumulated water deficit over the previous months (from 1 to 24). So what is returned is a list for each SPEI timescale (1 to 24 months).
- Computes stability metrics:
  - Resistance (ability to withstand drought/wet events). Computed as a deviation from the mean. This metric is computed as the formula defined in Isbell et al., Nature, 2015
  - Resilience (recovery after disturbance). Computed as the amount of biomass that recover one year after the drought. This metric is computed as the formula defined in Isbell et al., Nature, 2015
  - Log-ratio metrics
- Perform phylogenetic diversity analyses
- Perform trait analyses (PCA, trait representativity, CWM + FD along PCA axes)
- Merge everything together. *Note that we only provide for SPEI 3, 6, 9 in the git to reduce the size of the data. But running the analyses can generate the full list with the 24 SPEI timescales*

### Script 2: Analyse the data

Analyse the data along LUI, SPEI gradients, along each separate LUI component (fertilization, mowing, grazing) and then pairing resilience and resistance.
 

### Script 3: Making the figure

The figures are organized in three sections, based on the paper ones.
