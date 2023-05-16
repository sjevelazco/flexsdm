# flexsdm 1.3.3
- Improvements in `correct_colinvar` It is now possible to sample rasters to reduce machine memory and speed up the process
- Improvements in `sdm_predict` It is possible predict model in chunks to reduce machine memory
- `p_extra`, `p_pdp`, and `p_psp` were fixed 
- New function `p_psp` Partial Dependent Surface Plot
- New function `data_psp` Calculate data to construct partial dependence plots
- Improvements in `p_dpd` Calculate data to construct partial dependence plots

# flexsdm 1.3.2

- New function `p_extra` Graphical exploration of extrapolation or suitability pattern in the environmental and geographical space
- New function `p_pdp` Partial Dependent Plot
- New function `data_pdp` Calculate data to construct partial dependence plots
 

# flexsdm 1.3.1

- New argument "crs" added to function `msdm_posteriori`
- New argument "sp_name" in `sample_background` and `sample_pseudoabs`
- raster, flexclust, ape, and sp were removed from dependencies  
- Functions using CRS data have improved codes
- It is possible use a numeric value to specify threshold in `msdm_posteriori` 
- `extra_eval` can use tibble or SpatRaster object in env_calib argument 
- `extra_truncate` has a new argument to define values used for model truncation and documentation was improved. 
# 
