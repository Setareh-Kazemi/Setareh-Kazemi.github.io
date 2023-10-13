## Project Title

**Functional Regression for Upper Extremity Fatigue Analysis during Dynamic Order Picking**

--- 

## Project Members:  
- **Setareh Kazemi Kheiri, Department of Industrial and Systems Engineering, University at Buffalo**, **skazemik@buffalo.edu**
- **Sahand Hajifar, Department of Industrial and Systems Engineering, University at Buffalo**, **sahandha@buffalo.edu**
- **Linh Tran, Farmer School of Business, Miami University**, **tranlh2@miamioh.edu**
- **Zahra Vahedi, Department of Industrial and Systems Engineering, University at Buffalo**, **zahravah@buffalo.edu**
- **Hongyue Sun, Department of Industrial and Systems Engineering, University at Buffalo**, **hongyues@buffalo.edu**
- **Fadel M. Megahed, Farmer School of Business, Miami University**, **fmegahed@miamioh.edu**
- **Lora A. Cavuoto, Department of Industrial and Systems Engineering, University at Buffalo**, **loracavu@buffalo.edu**

---
## Objective:

In this document all the code, results and analysis related to the associated paper entitled "Explaining the Variability in the Profiles of Ratings of Perceived Exertion for a Dynamic Upper Extremity Task: A Functional Regression Approach" are provided. The main objective of doing this research was gaining a better understanding of how different factors contribute to the development of fatigue during manual material handling (MMH) operations in a simulated warehousing environment. Our approach is divided into the following main steps:

1.  Data pre-processing;

2.  Data transformation;

3.  Functional feature extraction;

4.  Functional regression on the transformed data with different set of scalar and functional features;

5.  Benchmark models;

6.  Checking models' assumptions (residual analysis)

![Image of Framework](functional%20regression-Phase-1.png)
![Image of Framework](functional%20regression-Phase-2.png)

---
## Introduction of the Experiment and Data:
In this repository you can find the code for the functional regression to explain the variability in the profiles of RPEs(Ratings of Perceived Exertion) which is a subjective fatigue indicator. Data used in this study is driven from an experiment conducted to assess the fatigue accumulation of **upper limb fatigue** in a repetitive overhead load lifting task. The experiment was conducted with four different conditions, with varying: (a) **task weights:** 1.5 and 2.5 kg, and (b) **task paces:** 5, 10, and 15 bpm. The four conditions are based on the following combinations of pace and weights: 

- 5 bpm -- 2.5 kg,   
- 10 bpm -- 2-5 kg,   
- 15 bpm -- 2.5 kg, and   
- 15 bpm -- 1.5 kg 

A total of 17 people participated in this experiment. Each session of the experiment consisted of three 45-minute periods, with two 15-minute breaks in between the work periods. The [Borg 0-10 Ratings of Perceived Exertion (RPE) Scale](https://my.clevelandclinic.org/health/articles/17450-rated-perceived-exertion-rpe-scale) were captured for each subject every 5 minutes of the experiment.  A demonstration of the different stages of the experiment task is shown in the below figure.

![Image of workstation](Simulated%20workstation.png)

---
## Data Files: 

In the [**data**](data) folder, there is an [**input**](input) folder which includes the below files:
   1. [Anthropometric_data.RData](data/inputs/Anthropometric_data.RData): In this file the anthropometrics of participants are saved including their: gender, age, height(cm), weight(kg), waist circumference (cm), hip circumference (cm), and body mass index (BMI).
   2. [Experimental_Design.xlsx](data/inputs/Experimental_Design.xlsx): In the sheet "For analysis" of this file the random order of task conditions (combinations of load and pace) are assigned to different participants. This file is later used in the analysis to understand what task condition was performed in each session of the experiment. 
   3. [Features.RData](data/inputs//Features.RData): This file contains the statistical feaures extracted from five-bottle cycles of movement in the experiment.
   4. [IMU_Raw.txt](data//inputs/IMU_Raw.txt): This file contains the link to access the large raw IMU data, provided in gogle drive for download.
   5. [corrected_Changepoints.csv](data//inputs/corrected_Changepoints.csv): In this file the changepoints identifying the start and end of the first 45-minute of experiment are provided.
   6. [y.RData](y.RData): In this file the raw RPE scores for the 44 full sessions are provided. full sessions are related to those sessions which the participants were able to finish the whole 45-minute duration of the experiment.
   7. [y_T_val.RData](y_T_val.RData): In this file the transformed RPE scores for the 44 full sessions are provided which is calculated using the multivariate BOX Cox method. The codes for extracting the RPE and TRPE values are provided in our previous project which can be found in ['functional_anova'](./functional_anova).
  

---
## Codes and Results: 

The codes related to all analyses represented in this study can be found in the ['functional_regression_final_SK.qmd'](functional_regression_final_SK.qmd) file which can be opened and rendered in R Studio after the data files are stored in a local directory same as this file. The rendered 'HTML' file can also be found in this repository entitled ['functional_regression_final_SK.html'](functional__regression_final_SK.html). The plots and figures related to different parts of this analysis can be found in the [**figs**](figs) folder.

