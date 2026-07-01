# Pupillometry-Based Data Analysis for Diagnostic Applications

This repository contains the code developed for the Master's thesis **"Pupillometry-Based Data Analysis for Diagnostic Applications"** by **Frederik N. H. Lundgren** and **Weinan Xiong**.

The project was completed as part of the M.Sc. programme at the **Technical University of Denmark (DTU)** in collaboration with **Rigshospitalet, Copenhagen University Hospital**.

## Project overview

The thesis investigates three interconnected aspects of quantitative pupillometry in patients with disorders of consciousness (DoC):

1. **NPi reconstruction**
   - Reconstruction of the proprietary Neurological Pupillary Index (NPi) using interpretable statistical and machine learning models.

2. **Consciousness modelling**
   - Investigation of the relationship between pupillary light reflex (PLR), light-off response (LOR), and clinical consciousness scores using both cross-sectional and longitudinal analyses.

3. **Survival prediction**
   - Evaluation of the prognostic value of pupillometric features and clinical scores for predicting 90-day survival.

## Repository structure

The repository is organized into folders corresponding to the three main research topics. Each folder contains the primary Jupyter notebooks used in the analyses.

Additional exploratory analyses, supplementary notebooks, and intermediate work can be found in the `Final/` directory.

## Data availability

The data used in this project were collected as part of an ongoing clinical study at Rigshospitalet.

Due to Danish and EU regulations regarding patient privacy and sensitive health information, the dataset **cannot be shared publicly**. Access to the data is therefore only available through Rigshospitalet and the associated research project.

For the same reason, many notebook outputs, print statements, and data visualizations have been removed or limited in this repository.

## Reproducing the analyses

The analyses were developed using the following environment:

```
Python version : 3.13.7
Operating system : Windows 11
OS version : 10.0.22631
Architecture : AMD64
Processor : Intel64 Family 6 Model 158 Stepping 10, GenuineIntel
```

All required Python packages are listed in `requirements.txt`.

## Notes

Because the original dataset is not publicly available, the notebooks cannot be executed end-to-end without access to the clinical data. Nevertheless, the repository documents the complete analysis pipeline, modelling approaches, and implementation used throughout the thesis.

## Authors

- Frederik N. H. Lundgren
- Weinan Xiong

Technical University of Denmark (DTU)

In collaboration with Rigshospitalet, Copenhagen University Hospital.
