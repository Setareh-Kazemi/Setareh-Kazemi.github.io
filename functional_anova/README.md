# Functional_anova_niosh
--- 

## Project Members:  
- **Setareh Kazemi Kheiri, Department of Industrial and Systems Engineering, University at Buffalo**, **skazemik@buffalo.edu**
- **Zahra Vahedi, Department of Industrial and Systems Engineering, University at Buffalo**, **zahravah@buffalo.edu**
- **Hongyue Sun, Department of Industrial and Systems Engineering, University at Buffalo**, **hongyues@buffalo.edu**
- **Fadel M. Megahed, Farmer School of Business, Miami University**, **fmegahed@miamioh.edu**
- **Lora A. Cavuoto, Department of Industrial and Systems Engineering, University at Buffalo**, **loracavu@buffalo.edu**

---
## Objective:

This document provides the code, results, and analysis for evaluating the development of fatigue during manual material handling (MMH) operations in a simulated warehousing environment. Our approach is divided into the following main steps: 

1. Data pre-processing;

2. Data transformation;

3. Functional ANOVA on the transformed data;

4. Data clustering; and

5. Functional ANOVA on the clustered transformed data

---
## Introduction of the Experiment and Data:
In this repository you can find the code for the functional ANOVA to find the significant factors impacting the objective and subjective fatigue indicator. Data used in this study is driven from an experiment conducted to assess the fatigue accumulation of **upper limb fatigue** in a repetitive overhead load lifting task. The experiment was conducted with four different conditions, with varying: (a) **task weights:** 1.5 and 2.5 kg, and (b) **task paces:** 5, 10, and 15 bpm. The four conditions are based on the following combinations of pace and weights: 

- 5 bpm -- 2.5 kg,   
- 10 bpm -- 2-5 kg,   
- 15 bpm -- 2.5 kg, and   
- 15 bpm -- 1.5 kg 

A total of 17 people participated in this experiment. Each session of the experiment consisted of three 45-minute periods, with two 15-minute breaks in between the work periods. The [Borg 0-10 Ratings of Perceieved Exertion (RPE) Scale](https://my.clevelandclinic.org/health/articles/17450-rated-perceived-exertion-rpe-scale) were captured for each subject every 5 minutes of the experiment. Furthermore, isometric strength tests were taken every 9 minutes. In this analysis, only the data related to the first 45 minutes of experiments were used to assess the impact of task characteristics on the accumulated fatigue, as well as comparing the two fatigue indicators (RPE and muscle strength).

---
## Data Files: 


