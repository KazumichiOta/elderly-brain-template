# README: Elderly Brain Template Data and Metrics

1. Overview

This repository contains supplementary data related to the Elderly Brain Template developed and evaluated in our study, “Developing and Validating an Elderly Brain Template: A Comprehensive Comparison with MNI152 for Age-Specific Neuroimaging Analyses.” Specifically, it includes:
- Subcortical Dice coefficients for 315 participants from the IXI dataset.
- Whole-brain metrics (3D-SSIM, CC, MSE) for the same participants.
- Demographics (ID, sex, and age) of older adults (65–84 years) selected from the OASIS-1 dataset, matching the 2020 Japanese census statistics.

2. File Contents
2.1 IXI_region_dice.csv
- Description: Subcortical Dice coefficients for each IXI subject, comparing registration to the Elderly template or MNI152.
- Columns:
  - ID: Subject ID from IXI
  - Age: Actual age in years
  - Age Group: Categorical label (A=20–29, B=30–39, C=40–49, D=50–59, E=60–69, F=≥70) used for visualization only
  - template: "Elderly" or "MNI152"
  - Region: e.g., amygdala, brainstem, hippocampus, pallidum, etc.
  - Dice: Dice coefficient for that region

2.2 IXI_metrics(SSIM_CC_MSE).csv
- Description: Whole-brain metrics (3D-SSIM, CC, MSE) for each IXI subject under Elderly vs. MNI152 registration.
- Columns:
  - ID, Age, Age Group, template
  - Metric: "3D-SSIM", "CC", or "MSE"
  - Value: Numerical value of the given metric

2.3 Oasis_Japan_demographics.csv
- Description: Demographic information (ID, sex, and age) for older adults (65–84 years) extracted from OASIS-1, reflecting age and sex distributions from the 2020 Japanese census.
- Columns:
  - ID: OASIS-1 subject ID
  - Sex: "F" or "M"
  - Age (years): Actual age

3. Usage Guidelines
3.1 Analysis Pipeline
- Please refer to our Supplementary Materials or main paper for details on the software (ANTs, FreeSurfer, etc.) used to compute these metrics, including the exact versions and parameters.

3.2 Age Group Note
- The Age Group column (A, B, C, D, E, F) is intended purely for visualization (e.g., color coding in figures). In the actual statistical analyses, age was treated as a continuous variable.

3.3 Reanalysis
- CSV files can be loaded into the analysis environment (Python, R, etc.) to replicate or extend our findings regarding template performance across different ages.

4. License

This dataset is released under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) license.
- Attribution: Please cite our original paper when using these data.
- NonCommercial: Commercial use of these data is not permitted.

5. Citation

If you use or reference these data, please cite the following paper accordingly:

Kazumichi Ota*, Yoshihiko Nakazato, Genko Oyama.  
“Developing and Validating an Elderly Brain Template: A Comprehensive Comparison with MNI152 for Age-Specific Neuroimaging Analyses.”  
2024. doi:xxx/yyy

*Note: This reference is provisional. If the paper has not yet been published or assigned a DOI, please update the citation once the final publication details become available.

6. Contact
- Corresponding Author: Kazumichi Ota*
- Email: Kota24@saitama-med.ac.jp
If you have any questions, feedback, or encounter issues with these files, please feel free to contact us.

7. Acknowledgments
- We thank the IXI and OASIS projects for providing open-access MRI data.
- No specific funding was received for this study.