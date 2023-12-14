# music2Param

`lute` param class, method, and generic definitions for the MuSiC2 deconvolution algorithm. See `?music2Param-class` for more information.

## Install

Install this package from GitHub using:

```
devtools::install_github("metamaden/music2Param")
```

## Dependencies

This param class requires the `MuSiC2` software to run (available from GitHub at [xuranw/MuSiC](https://github.com/xuranw/MuSiC) and ).

To run the `MuSiC2` implementation, you will need to have installed an older version of `MuSiC` (< v1.0.0) as of `MuSiC2` v0.1.0.

A YML file has been included to set up a conda environment to run the main dependencies `./inst/yml/music2.yml`.

## Citations

Fan, Jiaxin. MuSiC2: cell type deconvolution for multi-condition bulk RNA-seq data. (2023) GitHub, 
R package version 0.1.0. URL: https://github.com/Jiaxin-Fan/MuSiC2.
 
Wang, Xuran and Jiaxin Fan. MuSiC: Multi-subject single cell deconvolution. (2022) GitHub, R package 
version 1.0.0. URL: https://github.com/xuranw/MuSiC.
 
Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution with multi-subject single-cell 
expression reference. Nat Commun 10, 380 (2019). https://doi.org/10.1038/s41467-018-08023-x.

Jiaxin Fan, Yafei Lyu, Qihuang Zhang, Xuran Wang, Mingyao Li, Rui Xiao, MuSiC2: cell-type deconvolution 
for multi-condition bulk RNA-seq data, Briefings in Bioinformatics, Volume 23, Issue 6, November 2022, 
bbac430, https://doi-org.proxy1.library.jhu.edu/10.1093/bib/bbac430.
