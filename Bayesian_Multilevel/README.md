## Project Title

**A Hierarchical Bayesian Modeling Approach to Predict Upper Extremities' Fatigue in a Dynamic Order-Picking Task**

--- 

## Project Members:  
- **Setareh Kazemi Kheiri, Department of Industrial and Systems Engineering, University at Buffalo**, **skazemik@buffalo.edu**
- **Hongyue Sun, College of Engineering, University of Georgia**, **hongyuesun@uga.edu**
- **Fadel M. Megahed, Farmer School of Business, Miami University**, **fmegahed@miamioh.edu**
- **Lora A. Cavuoto, Department of Industrial and Systems Engineering, University at Buffalo**, **loracavu@buffalo.edu**

---
## Objective:

In this document all the code, and inputs related to the associated paper entitled "A Hierarchical Bayesian Modeling Approach to Predict Upper Extremities' Fatigue in a Dynamic Order-Picking Task" are provided. The main objective of doing this research was prediction of fatigue during manual material handling (MMH) operations in a simulated warehousing environment. Our approach is divided into the following main steps:

1.  Data pre-processing;

2.  Functional feature extraction;

3.  Fitting hierarchical regression to predict fatigue;

5.  Evaluation of models' prediction;


![Image of Framework](Bayesian%20Regression-Framework.png)


---
## Introduction of the Experiment and Data:
The experiment in this project is the same as the one explained in the [functional regression project](./functional_regression). 
The median features extracted from the wearable sensors are used as inputs to the models trained in this projects.

---
## Data Files: 

The [**data**](data) folder includes the below files:
   1. [Anthrpmtrc_data.RData](data/Anthrpomtrc_data.RData): In this file the anthropometrics of participants are saved including their: gender, age, height(cm), weight(kg), waist circumference (cm), hip circumference (cm), and body mass index (BMI).
   2. [input_data_sensors.RData](data/input_data_sensors.RData): This file contains the median feaures extracted from werable sensors attached to the trunk, wrist, and upperarm of participants.
   3. [y.RData](y.RData): In this file the raw RPE scores for the 43 full sessions are provided. full sessions are related to those sessions which the participants were able to finish the whole 45-minute duration of the experiment.
   4. [main_df.RData](main_df.RData): In this dataframe the information from task factors, individual characteristics, and sensors' median features are aggregated for all 43 participants. This file is needed to run the [metrics&plots_RPE.Rmd](metrics&plots_RPE.Rmd) file.
Additionally, the fitted models are stored in [**this Figshare link**](https://doi.org/10.6084/m9.figshare.28447265.v1) for easier access to run the [metrics&plots_RPE.Rmd](metrics&plots_RPE.Rmd). 
  
---
## Codes and Results: 

The codes related to the main analyses represented in this study can be found in the ['hierarchical_regression_brms.Rmd'](hierarchical_regression_brms.Rmd) file which can be opened and rendered in R Studio after the data files are stored in a local directory same as this file. The codes related to calculation of evaluation metrics and the prediction plots can be found in [metrics&plots_RPE.Rmd](metrics&plots_RPE.Rmd). The plots and figures related to different parts of this analysis can be found in the [**figs**](figs) folder.

