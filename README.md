## **LearningToLearn**

This repository contains MATLAB routines used to generate the Figures of Levi, Aviv, and Stark (2024) "Learning to learn: Single session acquisition of new rules by freely moving mice", PNAS Nexus (https://academic.oup.com/pnasnexus/advance-article/doi/10.1093/pnasnexus/pgae203/7676433)

## Overview
Humans excel at learning from examples and adapting to new rules, but conditions for
efficient learning in non-human subjects are unclear. We explored the rapid adaptability of mice
to new rules using a visual discrimination task. We found that mice can learn a new rule within a
single session, a capacity that enhances with experience and varies with rule difficulty.
Furthermore, mice exhibit flexibility in learning strategy, based on the physical conditions of the
task. 

The code available in this repository was used to perform the analysis resulting in Figures 2-4 and supplementary Figures S1-S2 of the paper.

## To run the code
1. Download all routines and data.
2. Open MATLAB and add the repository folder to the path.
3. Use the Matlab mex.m function:
   - mex parsec.c
   - mex zsc.c
4. Download the libsvm package:
Libsvm is a simple, easy-to-use, and efficient software for SVM
classification and regression. It solves C-SVM classification, nu-SVM
classification, one-class-SVM, epsilon-SVM regression, and nu-SVM
regression. It also provides an automatic model selection tool for
C-SVM classification.
Libsvm is available at
http://www.csie.ntu.edu.tw/~cjlin/libsvm
5. Follow the Libsvm README file for instruction on the package installation
6. Run the routines. For example, to plot Figure 2, write in MATLAB:
- tablename = 'PATH_TO_DEPOSITORY_FOLDER\levi2024_dataset';
- LtL_make_Figure2(tablename);
 
## Routines

### Wrappers
- LtL_make_Figure2.m
  - Generates Figure 2
- LtL_make_Figure3.m
  - Generates Figure 3
- LtL_make_Figure4.m
  - Generates Figure 4
- LtL_make_FigureS1.m
  - Generates Figure S1
- LtL_make_FigureS2.m
  - Generates Figure S2

### Utilities
- LtL_bayescount.m
  - Compute the number of bins for bias calculation, using the Bayes counting procedure of Panzeri&Treves 96
- LtL_calc_bias_PT.m
  - Calculate an analytical estimate of MI bias from a joint probability table
- LtL_mutinf.m
  - Mutual information from empirical distributions/counts
- myjet.m
  - Modified jet with extreme values in pure R,B.
- LtL_ParseArgPairs.m
  - Flexible argument assigning
- parse.m
  - Parse a vector to sets of points
- imagescnan.m
  - Image in myjet colormap on current plot
- LtL_make_equal_bins
  - assign indices to data such that equal sized bins
 
### C Utilities
These utility functions are operating system specific and require compiling in MATLAB (see section 'To run the code' above). 
- parsec.c
  - Parse a sorted array
- zcs.c
  - Compute the Z-score of the columns of a given matrix

### Data
- levi2024_dataset - excell sheet with the mice data
  - Sheet1: Data for Figures 1-3, S1-S2
    - Reward location
    - Trial result
    - Mouse number
    - Block ordinate within session
    - Trial ordinate within block
    - Rule number
    - Rule ordinate for each mouse
    - session ordinate for each mouse
    - Intra-criterion distance
    - Inter-criterion distance
    - Success rate in the first trial of the session
    - Success rate in the first two trials of the session
    - Success rate in the first three trials of the session
    - Success rate in the first four trials of the session
    - Success rate in the first five trials of the session
    - Success rate in the first six trials of the session
    - Success rate in the first seven trials of the session
    - Success rate in the first eight trials of the session
    - Success rate in the first nine trials of the session
    - Success rate in the first ten trials of the session
    - SSL index
    - MSL index
    - First session leanred index
  - Sheet2: Data for Figure 4
    - Reward location
    - Trial result
    - Mouse number
    - Block ordinate within session
    - Trial ordinate within block
    - Rule number
    - Rule ordinate for each mouse
    - session ordinate for each mouse
    - Intra-criterion distance
    - Inter-criterion distance
    - SSL index
    - MSL index



