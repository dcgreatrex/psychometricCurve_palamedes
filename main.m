function main()
%-----------------------------------------------------------------------------------------------------
% Author: David Greatrex, PhD Candidate, University of Cambridge.  
% Date: 09/03/2016
% Name: Fitting psychometric functions using the Palamedes Matlab toolbox
%-----------------------------------------------------------------------------------------------------
clear all

%--------------------------------------------
% User input variables
%--------------------------------------------
% number of participants in the dataset
pNo = 20;   
% path to the main downloaded psychometricCurve_palamedes folder
mainPath = '/Users/dcg/Code/Data_and_Analysis/PhD/Experiments/4_accumulation/code/palamedes_curvefitting';
% path to the main data folder
dataPath = '/Users/dcg/Code/Data_and_Analysis/PhD/Experiments/4_accumulation/data/data_for_R';
% binary indicator: simulate standard error statistics? (1 = Y, 0 = N) 
sim_std_error = 0;
% number of simulations to compute standard error statistic
std_err_reps = 500; 
% binary indicator: simulate goodness of fit statistics? (1 = Y, 0 = N) 
sim_goodness_fit = 0;
% number of simulations to compute goodness of fit statistics
good_fit_reps = 500; 

% column name of the response time column in the dataset
c_name_rt = 'rt';
% column name of the primary factor column in the dataset (i.e. for which a separate model will be fit to each level)
c_name_mainFactor = 'periodicity';
% column name of the level factor on which to regress the psychometric curves
c_name_level = 'level';

% binary indicator: plot fitted psychometric functions on data and save to plot folder.
plot_gen = 1;

f1  = 'pNo';                 v1 = pNo;
f2  = 'mainPath';            v2 = mainPath;
f3  = 'dataPath';            v3 = dataPath;
f4  = 'sim_std_error';       v4 = sim_std_error;
f5  = 'std_err_reps';        v5 = std_err_reps;
f6  = 'sim_goodness_fit';    v6 = sim_goodness_fit;
f7  = 'good_fit_reps';       v7 = good_fit_reps;
f8  = 'c_name_rt';           v8 = c_name_rt;
f9  = 'c_name_mainFactor';   v9 = c_name_mainFactor;
f10 = 'c_name_level';        v10 = c_name_level;
f11 = 'plot_gen';            v11 = plot_gen;

meta_data = struct(f1, v1, f2, v2, f3, v3, f4, v4, ...
                   f5, v5, f6, v6, f7, v7, f8, v8, ...
                   f9, v9, f10, v10, f11, v11);
    

%--------------------------------------------
% Fir psychometric curves
%--------------------------------------------
fitting(meta_data)