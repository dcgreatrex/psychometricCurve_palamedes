# Fitting psychometric functions using the Palamedes toolbox
Author: David Greatrex, PhD Candidate, University of Cambridge.  
Date: 01/03/2016 -- Language: Matlab. -- Modifications:

This script uses functions taken from the Palamedes Matlab toolbox:
* Prins, N. & Kingdom, F.A.A. Palamedes: Matlab routines for analyzing psychophysical data. www.palamedestoolbox.org

## Summary:
This package provides an example of how to fit the cumulative distribution function to
psychoacoustic data collected during psychological experimentation at the University of
Cambridge.

## Data:
The data is taken from an experiment investigating the effects of rhythmic expectation on
complex decision-making. This used the psychophysical method of constant stimuli, in which
the target stimulus was varied systematically around a detection threshold.

## Purpose:
The purpose of the script is to fit the cumulative distribution function to each participant's 
data, individually by each condtion of interest. Slope and threshold values are then extracted 
from each fit and output to a table. Parametric bootstrapping is then applied to the 
fitted function as means of computing standard error values for both slope and threshold as well 
as computing goodness of fit statistics. These additions are found in the output table.

## Method:
The psychometric data from each participant and condition were fitted with sigmoidal cumulative 
distribution functions, each defined by three parameters: threshold, slope and lapse-rate as 
implemented in Palamedes toolbox (Prins and Kingdom, 2009). This is the same method used by
Rohenkohl et al, 2012 & Cravo et al, 2013. Guess rates were fixed at 0 (i.e., the participant 
always pressed the same response key) across subjects and conditions, and the three parameters 
were fit separately for each subject and condition.
* http://www.palamedestoolbox.org/overview.html
* http://www.palamedestoolbox.org/understandingfitting.html

## To run:
1. Download the content of the folder locally on to your hardrive including the folder structure. 
2. Open the Matlab script 'main.m'. 
3. Change the path on line 7 to represent the folder location in which you saved the downloaded content.
4. Run 'main.m'. All plots will then be saved into the downloaded filestructure plot folder.
Note - due to bootstrapping the script will take a long time to complete. ~20 minutes. You can check progress
by looking at how many plots have been output into the plot folder. There are 17 participants in total.

## References:
* [Prins, N. & Kingdom, F.A.A. Palamedes: Matlab routines for analyzing psychophysical data.](http://www.palamedestoolbox.org/index.html)

* [Rohenkohl, G. Cravo, A. Wyart, V. and Nobre, A. (2012). Temporal Expectation Improves the Quality of Sensory Information. The Journal of Neuroscience, 32(24), 8424–8428.](https://www.jneurosci.org/content/32/24/8424.full)

* [Cravo, A. Rohenkohl, G. Wyart, V. and Nobre, A. (2013). Temporal Expectation Enhances Contrast Sensitivity by Phase Entrainment of Low-Frequency Oscillations in Visual Cortex. The Journal of Neuroscience, 33 (9), 4002–4010.](https://www.jneurosci.org/content/33/9/4002.full)