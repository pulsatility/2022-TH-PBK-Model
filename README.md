# A Minimal Human Physiologically Based Kinetic (PBK) Model of Thyroid Hormones and Chemical Disruption of Plasma Thyroid Hormone Binding Proteins
Anish D. Bagga<sup>1</sup>, Brian P. Johnson<sup>2</sup>, and Qiang Zhang<sup>3</sup>

Published in  Frontiers in Endocrinology (2023) 14:1168663 (https://doi.org/10.3389/fendo.2023.1168663)

1. Emory College of Arts and Sciences, Emory University, Atlanta, GA 30322, USA

2. Department of Pharmacology and Toxicology, Michigan State University, East Lansing, MI 48824, USA

3. Gangarosa Department of Environmental Health, Rollins School of Public Health, Emory University, GA 30322, USA

  
**Abstract:**
The thyroid hormones (THs), thyroxine (T4) and triiodothyronine (T3), are under homeostatic control by the hypothalamic-pituitary-thyroid axis and plasma TH binding proteins (THBPs), including thyroxine-binding globulin (TBG), transthyretin (TTR), and albumin (ALB).  THBPs buffer free THs against transient perturbations and distribute THs to tissues. TH binding to THBPs can be perturbed by structurally similar endocrine-disrupting chemicals (EDCs), yet their impact on circulating THs and health risks are unclear. In the present study, we constructed a human physiologically based kinetic (PBK) model of THs and explored the potential effects of THBP-binding EDCs. The model describes the production, distribution, and metabolism of T4 and T3 in the Body Blood, Thyroid, Liver, and Rest-of-Body (RB) compartments, with explicit consideration of the reversible binding between plasma THs and THBPs. Rigorously parameterized based on literature data, the model recapitulates key quantitative TH kinetic characteristics, including free, THBP-bound, and total T4 and T3 concentrations, TH productions, distributions, metabolisms, clearance, and half-lives. Moreover, the model produces several novel findings. (1) The blood-tissue TH exchanges are fast and nearly at equilibrium especially for T4, providing intrinsic robustness against local metabolic perturbations. (2) Tissue influx is limiting for transient tissue uptake of THs when THBPs are present. (3) Continuous exposure to THBP-binding EDCs does not alter the steady-state levels of THs, while intermittent daily exposure to rapidly metabolized TBG-binding EDCs can cause much greater disruptions to plasma and tissue THs. In summary, the PBK model provides novel insights into TH kinetics and the homeostatic roles of THBPs against thyroid disrupting chemicals.

**Keywords:** Thyroid hormone, Thyroid hormone binding protein, Thyroxine-binding globulin, Transthyretin, Albumin, PBK model, Endocrine disrupting chemicals


#  MATLAB Code
- TH_PBK_Nonspatial_CMD.m: Main MATLAB code to generate steady-state results.
- TH_PBK_Nonspatial_CMD_tracer.m: MATLAB code to generate TH tracer simulation-related results.
- TH_PBK_Nonspatial_ODE.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD.m.
- TH_PBK_Nonspatial_ODE_clamped.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD.m for sensitivity analysis.
- TH_PBK_Nonspatial_ODE_tracer.m: MATLAB ODE code to be called by TH_PBK_Nonspatial_CMD_tracer.m for TH tracer simulations.
