#  MATLAB Code for a Minimal Human PBK Model of Thyroid Hormones (nonspatial)
Published in Bagga M., Johnson B., and Zhang Q. (2023) A Minimal Human Physiologically Based Kinetic (PBK) Model of Thyroid Hormones and Chemical Disruption of Plasma Thyroid Hormone Binding Proteins (accepted)

Purpose of files:
- TH_PBK_Nonspatial_CMD.m: Main MATLAB code to generate steady-state results.
- TH_PBK_Nonspatial_CMD_tracer.m: MATLAB code to generate TH tracer simulation-related results.
- TH_PBK_Nonspatial_ODE.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD.m.
- TH_PBK_Nonspatial_ODE_clamped.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD.m for sensitivity analysis.
- TH_PBK_Nonspatial_ODE_tracer.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD_tracer.m for TH tracer simulations.
