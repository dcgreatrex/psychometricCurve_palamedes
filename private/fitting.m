function fitting(meta_data)

%--------------------------------------------
para.thresh = [];
para.slope = [];
stderror.thresh = [];
stderror.slope = [];
gof.deviance = [];
%--------------------------------------------
% Output dataset structure
%--------------------------------------------
VarNames = {'subNo','Threshold', 'Slope', 'T.SdErr', 'S.SdErr','Deviance', 'goodFitPvalue'};
periodicResults = dataset([],[],[],[],[],[],[],'VarNames',VarNames);
aperiodicResults = dataset([],[],[],[],[],[],[],'VarNames',VarNames);
for k = 1:meta_data.pNo
    
    %--------------------------------------------
    % Load data
    %--------------------------------------------
    file = strcat('/',num2str(k),'/',num2str(k),'_edit.txt');
    wd = strcat(meta_data.dataPath,file);
    p = importdata(wd);
    colheaders = strsplit(p.textdata{1,:});
    colheaders = strrep( colheaders,'"','');
    c =mat2dataset(p.data, 'VarNames', colheaders);
    
    %--------------------------------------------
    % Summary statistics
    %--------------------------------------------
    rt = double(c(:, meta_data.c_name_rt));
    mu = mean(rt);
    sd = std(rt);

    %--------------------------------------------
    % Data cleaning
    %--------------------------------------------
    %upper_limit = mu + sd; % remove any trials in which the RT was greater than +3 SD from the mean.
    %c2 = c(rt < upper_limit,:);
    c2 = c 

    %--------------------------------------------
    % Subset data by experimental condition
    %--------------------------------------------
    main_factor = double(c2(:, meta_data.c_name_mainFactor))
    condition_1 = c2(main_factor == 1,:);
    condition_2 = c2(main_factor == 0,:);

    %--------------------------------------------
    % Input variables to the curve fitting loop
    %--------------------------------------------
    dsArray = {condition_1, condition_2};
    ProportionCorrectObserved = {};
    StimLevelsFineGrain = {};
    ProportionCorrectModel = {};
    parameters = {};
    standardErrors = {};
    goodnessFit = {};
    goodnessFitP = {};
    
    %--------------------------------------------
    % Curve fitting loop for both experimental conditions
    %--------------------------------------------
    for j = 1:length(dsArray)

        % define data
        inputData = dsArray{j};

        % Stimulus intensities
        level = double(inputData(:, meta_data.c_name_level));
        StimLevels = unique(level)';

        % Aggregated data
        NumPos = [];
        NumPosRight = [];
        OutOfNum = [];
        for i = 1:length(StimLevels)
            % Number of left responses (e.g., 'answer == -1' at each of the 
            % entries of 'StimLevels'
            tmp = inputData(level == StimLevels(i),:);
            out = sum(tmp.answer(:) == 1);
            NumPos = [NumPos, out];
            outLeft = sum(tmp.answer(:) == -1);
            NumPosRight = [NumPos, outLeft];
            % Number of trials at each entry of 'StimLevels'
            % OutOfNum = [OutOfNum, sum(inputData.testLevel(:) == StimLevels(i))];
            OutOfNum = [OutOfNum, sum(out + outLeft)];
        end
        disp(StimLevels);
        disp(NumPos);
        disp(OutOfNum);
        
        %Use the cumulative normal function following Knoblauch et al. (2012). Modelling psychophysical data in R. Springer
        PF = @PAL_CumulativeNormal   % PAL_CumulativeNormal   %PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull,
                             %PAL_Quick, PAL_logQuick, %PAL_CumulativeNormal, PAL_HyperbolicSecant

        %Threshold and Slope are free parameters, guess and lapse rate are fixed
        paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter

        %Parameter grid defining parameter space through which to perform a
        %brute-force search for values to be used as initial guesses in iterative
        %parameter search.
        searchGrid.alpha = min(StimLevels):.001:max(StimLevels);
        searchGrid.beta = 10.^[-1:.025:2];  %logspace(-1,3,101);
        searchGrid.gamma = 0;  %scalar here (since fixed) but may be vector
        searchGrid.lambda = 0:.01:0.08;
        lapseLimits = [0 1];

        %Perform fit
        disp('Fitting function.....');

        %Case of field names must match exactly. E.g., options.tolx = 1 e-12 will 
        %be ignored (no warning will be issued). Options structure must be passed
        %to PAL_PFML_Fit (see call to PAL_PFML_Fit below) or will go ignored.
        options.TolX = 1e-12;   %increase desired precision (by lowering tolerance)
        options.TolFun = 1e-12;
        options.MaxIter = 2000;
        options.MaxFunEvals = 2000;
        [paramsValues LL exitflag] = PAL_PFML_Fit(StimLevels,NumPos, ...
            OutOfNum,searchGrid,paramsFree,PF,'lapselimits',lapseLimits, 'searchOptions',options);
        disp('done:')
        message = sprintf('Threshold estimate: %6.4f',paramsValues(1));
        disp(message);
        message = sprintf('Slope estimate: %6.4f\r',paramsValues(2));
        disp(message);
        
        parameters{j} = paramsValues;
        
        %--------------------------------------------
        % Simulate standard error statistics
        %--------------------------------------------
        if(meta_data.sim_std_error)
            %Number of simulations to perform to determine standard error
            B=20;                  
            ParOrNonPar = 1;
            disp('Determining standard errors.....');
            if ParOrNonPar == 1
                [SD paramsSim LLSim converged] = PAL_PFML_BootstrapParametric(...
                    StimLevels, OutOfNum, paramsValues, paramsFree, B, PF, ...
                    'searchGrid', searchGrid);
            else
                [SD paramsSim LLSim converged] = PAL_PFML_BootstrapNonParametric(...
                    StimLevels, NumPos, OutOfNum, [], paramsFree, B, PF,...
                    'searchGrid',searchGrid);
            end
            disp('done:');
            message = sprintf('Standard error of Threshold: %6.4f',SD(1));
            disp(message);
            message = sprintf('Standard error of Slope: %6.4f\r',SD(2));
            disp(message);
            standardErrors{j} = SD;
        else
            standardErrors{j} = [1 1];
        end
        
        %--------------------------------------------
        % Simulate goodness of fit statistics
        %--------------------------------------------
        if(meta_data.sim_goodness_fit)
        % Number of simulations to perform to determine Goodness-of-Fit
            B=20;
            disp('Determining Goodness-of-fit.....');
            [Dev, pDev, devSim] = PAL_PFML_GoodnessOfFit(StimLevels, NumPos, OutOfNum, ...
                paramsValues, paramsFree, B, PF, 'searchGrid', searchGrid);

            disp('done:');
            %Put summary of results on screen
            message = sprintf('Deviance: %6.4f',Dev);
            disp(message);
            message = sprintf('p-value: %6.4f',pDev);
            disp(message);
            devSim
            
            goodnessFit{j} = Dev;
            goodnessFitP{j} = pDev;
        else
            goodnessFit{j} = 1;
            goodnessFitP{j} = 1;
        end
        
        % Save output data to lists
        ProportionCorrectObserved{j} = NumPos./OutOfNum; 
        StimLevelsFineGrain{j} = [min(StimLevels):max(StimLevels)./1000:max(StimLevels)];
        ProportionCorrectModel{j} = PF(paramsValues,StimLevelsFineGrain{j}); 
    end
    
    %--------------------------------------------
    % Sompile simulated statistics for both of the periodicity conditions
    %--------------------------------------------
    pResultsTmp = dataset(k,parameters{1}(1), parameters{1}(2), standardErrors{1}(1), standardErrors{1}(2), ...
               goodnessFit{1}(1), goodnessFitP{1}(1), ...
              'Varnames', VarNames);
    apResultsTmp = dataset(k,parameters{2}(1), parameters{2}(2), standardErrors{2}(1), standardErrors{2}(2), ...
               goodnessFit{2}(1), goodnessFitP{2}(1), ...              
               'Varnames', VarNames);
     
    % update master output datasets      
    periodicResults = vertcat(periodicResults, pResultsTmp);
    aperiodicResults = vertcat(aperiodicResults, apResultsTmp);


    %--------------------------------------------
    % Plot fitted data for each participant
    %--------------------------------------------
    if(meta_data.plot_gen)
        
        for (i = 1:2)
            % plot data for each participant
            fig = figure('name', strcat('P_',num2str(k),' Maximum Likelihood Psychometric Function Fitting'));
            axes
            col = {[0 .7 0], [.7 0 0]};
            pLabels = {'Periodic', 'Aperiodic'};
            tsLabels = {'Threshold', 'Slope'};
            hold on   
            % add points
            x = StimLevels;
            y1point = ProportionCorrectObserved{1};
            y2point = ProportionCorrectObserved{2};
            plot(x, y1point, 'k.','color', col{1}, 'markersize', 40);
            plot(x, y2point, 'k.','color', col{2}, 'markersize', 40);
            % add legend
            legend('Periodic','Apeirodic','Location','northwest');
            hold on
            % add reference line at y = 0.5
            hline = refline([0 0.5]);
            set(hline,'LineStyle','--', 'color', [0.5 0.5 0.5]);
            hold on
            % add threshold lines
            if(i == 2)
                line([parameters{1}(1) parameters{1}(1)], [0,0.5], 'color', col{1}, 'linewidth', 2);
                line([parameters{2}(1) parameters{2}(1)], [0,0.5], 'color', col{2}, 'linewidth', 2);
                hold on
            end
            % add curve predicitons
            plot(StimLevelsFineGrain{1},ProportionCorrectModel{1}, '-', 'color', col{1}, 'linewidth', 4);
            plot(StimLevelsFineGrain{2},ProportionCorrectModel{2}, '-', 'color', col{2}, 'linewidth', 4);
            hold off

            % add threshold and slope values
            if(i == 2)
                m1 = strcat('Periodic  : T = ', num2str(parameters{1}(1)), ', S = ', num2str(parameters{1}(2)));
                m2 = strcat('Aperiodic: T = ', num2str(parameters{2}(1)), ', S = ', num2str(parameters{2}(2)));
                text(0, 0.2, m1, 'fontsize', 14);
                text(0, 0.1, m2, 'fontsize', 14);
            end
            % label axes
            set(gca, 'fontsize',16);
            set(gca, 'Xtick',StimLevels);
            axis([min(StimLevels) max(StimLevels) 0 1]);
            title(strcat(num2str(k),': Maximum Likelihood Fit'));
            xlabel('Displacement of probe tone from average');
            ylabel('% left responses');
            
            % save plot to file
            fileName = [meta_data.mainPath, '/plots/', num2str(k), '_plot ', num2str(i), '.jpeg'];
            saveas(fig,fileName);
        end
    end
end

%--------------------------------------------
% Save output datasets to file for analysis
%--------------------------------------------
outP = strcat(meta_data.mainPath, '/outTables/PeriodicFitValues.txt');
export(periodicResults,'file', outP)
outAP = strcat(meta_data.mainPath, '/outTables/AperiodicFitValues.txt');
export(aperiodicResults,'file', outAP)