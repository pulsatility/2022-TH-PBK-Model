#  MATLAB Code for Human PBK Models of Thyroid Hormones

(A) Purpose of files for nonspatial model published in Bagga M., Johnson B., and Zhang Q. (2023) A Minimal Human Physiologically Based Kinetic (PBK) Model of Thyroid Hormones and Chemical Disruption of Plasma Thyroid Hormone Binding Proteins. Frontiers in Endocrinology 14:1168663 (https://doi.org/10.3389/fendo.2023.1168663)
- TH_PBK_Nonspatial_CMD.m: Main MATLAB code to generate steady-state results.
- TH_PBK_Nonspatial_CMD_tracer.m: MATLAB code to generate TH tracer simulation-related results.
- TH_PBK_Nonspatial_ODE.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD.m.
- TH_PBK_Nonspatial_ODE_clamped.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD.m for sensitivity analysis.
- TH_PBK_Nonspatial_ODE_tracer.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD_tracer.m for TH tracer simulations.

(B) Purpose of files for spatial model:
- TH_PBK_Spatial_CMD.m: Main MATLAB code to generate steady-state results.
- TH_PBK_Spatial_ODE.m: MATLAB ODE code to be called by TH_PBK_Spatial_CMD.m.
- TH_PBK_Spatial_CMD_gradient.m: MATLAB code to generate steady-state results for parameters with gradients.
- TH_PBK_Spatial_ODE.m: MATLAB ODE code to be called by TH_PBK_Spatial_CMD_gradient.m.
- TH_PBK_Spatial_code_for_Fig_S1.m: MATLAB code to generate Fig. S1.
