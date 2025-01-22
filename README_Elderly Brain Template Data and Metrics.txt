README_Elderly Brain Template Data and Metrics

1. Overview

This repository contains supplementary data related to the Elderly Brain Template developed and evaluated in our study,
“Developing and Validating an Elderly Brain Template: A Comprehensive Comparison with MNI152 for Age-Specific Neuroimaging Analyses.”

Specifically, it includes:
	•	Subcortical Dice coefficients for 315 participants from the IXI dataset.
	•	Whole-brain metrics (3D-SSIM, CC, MSE) for the same participants.
	•	Demographics (ID, sex, and age) of older adults (65–84 years) selected from the OASIS-1 dataset, reflecting the 2020 Japanese census distribution.
	•	Iterative convergence metrics (MSE, CC, 3D-SSIM) capturing how the Elderly Template was refined over 20 iterative updates (phases).

2. File Contents

2.1 IXI_region_dice.csv
	•	Description: Subcortical Dice coefficients for each IXI subject, comparing registration to the Elderly template or MNI152.
	•	Columns:
	•	ID: Subject ID from IXI
	•	Age: Actual age in years
	•	Age Group: Categorical label (A=20–29, B=30–39, C=40–49, D=50–59, E=60–69, F=≥70) for visualization only
	•	template: "Elderly" or "MNI152"
	•	Region: e.g., amygdala, brainstem, hippocampus, pallidum, etc.
	•	Dice: Dice coefficient measuring subcortical overlap between the registered subject image and the chosen template

2.2 IXI_metrics(SSIM_CC_MSE).csv
	•	Description: Whole-brain metrics (3D-SSIM, CC, MSE) for each IXI subject under Elderly vs. MNI152 registration.
	•	Columns:
	•	ID, Age, Age Group, template
	•	Metric: One of "3D-SSIM", "CC", or "MSE"
	•	Value: Numerical value of that metric

2.3 Oasis_Japan_demographics.csv
	•	Description: Demographic information (ID, sex, age) for older adults (65–84 years) extracted from OASIS-1, matching the age and sex distributions reported in the 2020 Japanese census.
	•	Columns:
	•	ID: OASIS-1 subject ID
	•	Sex: "F" or "M"
	•	Age (years): Actual age

2.4 ElderlyTemplate_construction(SSIM_CC_MSE).csv
	•	Description:
This file contains the iterative metrics used to evaluate the convergence of our Elderly template construction. We performed a total of 20 updates (phases), and after every two steps, we checked 3D-SSIM (three-dimensional structural similarity index), CC (cross-correlation), and MSE (mean squared error) to determine whether the template had sufficiently converged.
	•	Columns:
	1.	Phase–Iteration: Labels like phase1-1, phase1-2, phase2-1, etc., indicating the phase and step number.
	2.	FileName: The warped subject image filename at that iteration (e.g., phase1_input0002-modality0-WarpedToTemplate.nii.gz).
	3.	CC (Cross-Correlation): A voxelwise measure of intensity correlation. Higher CC implies better alignment.
	4.	MSE (Mean Squared Error): The mean of squared intensity differences. Lower MSE indicates fewer discrepancies.
	5.	SSIM (3D-SSIM): Structural similarity index focusing on spatial patterns. Higher SSIM signals better structural match.

These iterative metrics reveal how each subject’s registration to the evolving Elderly template improved (or stabilized) at each step. In the main paper, Figure 2 provides a box-plot visualization of these changes in CC, MSE, and SSIM across the 20 updates.

3. Usage Guidelines

3.1 Analysis Pipeline

Please refer to our Supplementary Materials or the main paper for details on the software (ANTs, FreeSurfer, etc.) used to compute these metrics, including exact versions and parameters.

3.2 Age Group Note

The Age Group column (A, B, C, D, E, F) is intended purely for visualization (e.g., color-coding in figures). In the actual statistical analyses, age was treated as a continuous variable.

3.3 Reanalysis

By loading these CSV files into your preferred environment (Python, R, etc.), you can replicate or expand our findings on template performance across different ages—both in terms of subcortical overlap (Dice) and whole-brain metrics (3D-SSIM, CC, MSE). For the Elderly template construction process specifically, ElderlyTemplate_construction(SSIM_CC_MSE).csv documents each iteration’s alignment quality.

4. License

This dataset is released under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) license.
	•	Attribution: Please cite our original paper when using these data.
	•	NonCommercial: Commercial use is not permitted.

5. Citation

If you use or reference these data, please cite the following paper accordingly:

Kazumichi Ota*, Yoshihiko Nakazato, Genko Oyama.
“Developing and Validating an Elderly Brain Template: A Comprehensive Comparison with MNI152 for Age-Specific Neuroimaging Analyses.”
2024. doi:xxx/yyy

	Note: This reference is provisional. Once the paper is published or assigned a DOI, please update this citation with the final publication details.

6. Contact
	•	Corresponding Author: Kazumichi Ota*
	•	Email: Kota24@saitama-med.ac.jp

If you have any questions, feedback, or encounter issues with these files, please feel free to contact us.

7. Acknowledgments
	•	We thank the IXI and OASIS projects for providing open-access MRI data.
	•	No specific funding was received for this study.