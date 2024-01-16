%plotOptimiserResults.m This takes workspaces saved in the Data folder
%which contain the results of optimiser rsults and analyses and plots the
%data.

%Important variables to checkbefore running: Boundary, rows, cols, 
%Stochastic, CoupledCells; do they match the data being used from the 
%optimsier?

clear; clc; close all
addpath('Functions')
addpath('Data')

%% Choose the colour of all graphs
colors = magma(100);       %Plot colours
color1 = [0 0 0 0.05];
color2 = [colors(65,:) 0.05];

GraphAppearance = 0;       %0 = Normal, 1 = Dark

if GraphAppearance==1   
    BackgroundColour=38/255*[1 1 1]; TextColour=[1 1 1];
    get(0,'Factory');                           set(0,'defaultfigurecolor',BackgroundColour)
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor',BackgroundColour)
    set(0,'defaultAxesXColor',TextColour);      set(0,'defaultAxesYColor',TextColour)
    set(0,'defaultLegendTextColor',TextColour); set(0,'defaultTextColor',TextColour)  
    set(0,'defaultAxesZColor',TextColour)
elseif GraphAppearance==0
    get(0,'Factory');                           set(0,'defaultfigurecolor',[0.94 0.94 0.94])
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor','w')
    set(0,'defaultAxesXColor','k');             set(0,'defaultAxesYColor','k')
    set(0,'defaultLegendTextColor','k');        set(0,'defaultTextColor','k')
end

%% Define grid size
rows=10;                        %Number of rows of cells in hexagonal grid
cols=6;                         %Number of columns of cells in hexagonal grid

%% Define simulation time and step size
t0=0;                           % Start time
tf_hours=100;                   % Final time (hours)

%% Set up simulation inputs
%Important options!
Stochastic = 1;                 %1 = stochastic, 0 = deterministic
CoupledCells = 1;               %Simulate with Notch-Delta coupling = 1, uncoupled = 0
Boundary = 1;                   %0 = 'hard' boundary, 1 = periodic boundary (set to 1)
TurnOffAutorepression = 0;      %Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell = 1, with autorepression = 0

%% Load parameter sets identified by the optimiser
% load("Model 1 Uncoupled - Rows=10  Cols=6  Stochastic=1  Number of runs = 6000.mat") % Model 1: Uncoupled multicellular model - remember to set the 'CoupledCells' variable above to 0 if using these parameter sets
load("Model 2 Coupled - Rows=10  Cols=6  Stochastic=1  Number of runs = 6000.mat")     % Model 2: Coupled parameter set - remember to set the 'CoupledCells' variable above to 1 if using these parameter sets
N = numel(FVAL); %Number of identified parameter sets

%% Filter data
fval_cutoff      = -0.25;       %Filter for values less than this (nominal = -0.25)
coherence_cutoff = 1;           %Filter for values less than this (nominal = 1)
abundance_cutoff = 2000;        %Filter for values greater than this (nominal = 3000)
elemKeep=find(FVAL < fval_cutoff &  max(COHERENCE_MEAN, [], 2) < coherence_cutoff    &   ABUNDANCE_MEDIAN(:,1) > abundance_cutoff == true);

XMIN=XMIN(elemKeep, :);
FVAL=FVAL(elemKeep);
COHERENCE_MEAN =  COHERENCE_MEAN(elemKeep, :);
PERIOD_MEAN = PERIOD_MEAN(elemKeep, :);
COV = COV(elemKeep, :);

for elem = 1 : numel(elemKeep)
    COHERENCE_IND_CELLS_Filtered{elem, 1} = COHERENCE_IND_CELLS{elemKeep(elem)};
    PERIOD_IND_CELLS_Filtered{elem, 1} = PERIOD_IND_CELLS{elemKeep(elem)};
end


%% Order parameter sets by a certain output value
[FVAL,sortIdx] = sort(FVAL); %Reorder by largest FVAL 

XMIN = XMIN(sortIdx,:);
FVAL = FVAL(sortIdx);
COHERENCE_MEAN = COHERENCE_MEAN(sortIdx, :);
PERIOD_MEAN = PERIOD_MEAN(sortIdx, :);
COV = COV(sortIdx, :);

for idx = 1 : numel(sortIdx)
    COHERENCE_IND_CELLS_Filtered{idx, 1} = COHERENCE_IND_CELLS{sortIdx(idx)};
    PERIOD_IND_CELLS_Filtered{idx, 1} = PERIOD_IND_CELLS{sortIdx(idx)};
end

NumberOfParamSets = length(FVAL);
figNum = 10;
subPlotNum = 1;
modelParams={};


%% Calc median periods
for paramSet = 1:NumberOfParamSets
    PERIOD_MEDIAN(paramSet, :) = [median(PERIOD_IND_CELLS_Filtered{paramSet}(:,1)), median(PERIOD_IND_CELLS_Filtered{paramSet}(:,2))];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(201)
set(gcf,'renderer','painters') %For EPS file export
clf

subplot(132)
violinplot(PERIOD_MEDIAN,{'Normal protein deg', '1.1X degradation'}, 'QuartileStyle','shadow', 'ShowMean', true, 'ViolinColor', [color1(1:3); color2(1:3)]);
title('Average temporal period')
ylabel('Period (h)')
axis square
drawnow
set(gca, 'fontsize', 15)

subplot(133)
violinplot(COHERENCE_MEAN, {'Normal protein deg', '1.1X degradation'}, 'QuartileStyle','shadow', 'ShowMean', true, 'ViolinColor', [color1(1:3); color2(1:3)]); hold on
title('Average individual cell coherence')
axis square
set(gca, 'fontsize', 15)


%% Statistical tests
[pVal_Coherence, ~] = kruskalwallis(COHERENCE_MEAN);
[pVal_Period, ~] = kruskalwallis(PERIOD_MEDIAN);
% [pVal_CoV_indCells, ~] = kruskalwallis(CoV_avgIndCells);


%% Histogram plots
%Select values that have error below a certain value to plot
XMINplot=XMIN;
FVALplot=FVAL;

selectIdx = XMINplot(:,10)<100;
XMINplot=XMINplot(selectIdx, :);
FVALplot=FVALplot(selectIdx);

fprintf(sprintf('\n Percentage that are below error threshold = %.f', 100*length(elemKeep)/N))

axisLabels={'a_m', 'a_p', 'P_{H0}', 'n_H', 'Tau_H', 'P_{ND0}', 'n_{ND}', 'Tau{ND}', 'mRNA T_{1/2} ', 'protein T_{1/2}'};
figure(20)
clf
for m=1:10
    subplot(2,5,m)
%     histogram(XMINplot(:,m), 'FaceColor', 1/255*[156 213 153])
    histogram(XMINplot(:,m), 'FaceColor', 0.5*[1 1 1])
    xlim([LB(m) UB(m)])
    axis square
    xlabel(axisLabels(m))
    set(gca, 'FontSize', 12)
    
end

sgtitle('Parameter value distributions')
assignin('base','XMIN',XMIN)
assignin('base','FVAL',FVAL)

%Scatter plot
figure(3)
param1 = 10;
param2 = 3;
param3 = 6;
scatter3(XMINplot(:, param1), XMINplot(:, param2), XMINplot(:, param3), 40, -FVALplot, 'filled')
xlabel('Protein half-life (min)')
ylabel('P_{H0}')
zlabel('P_{ND0}')

colormap(inferno(255))
colorbar

for paramSet = 1:75:NumberOfParamSets

    clc
    fprintf(sprintf('Progress: %.f/%.f', paramSet, NumberOfParamSets));

    X = XMIN(paramSet,:);

    a_m   = X(1);         %mRNA production rate
    a_p   = X(2);         %Protein production rate
    P_H0  = X(3);         %Autorepression repression threshold
    n_H   = X(4);         %Autorepression Hill coefficient
    TauH  = X(5);         %Autorepression time delay
    
    P_ND0 = X(6);         %Lateral inhibition repression threshold
    n_ND  = X(7);         %Lateral inhibition Hil coefficient
    TauND = X(8);  
    
    u_m   = log(2)/X(9);  %mRNA degredation rate
    u_p   = log(2)/X(10); %Protein degredation rate
    
    
    %% Parameter set-up
    cells = rows*cols; %Total number of cells
    tf = tf_hours*60;  % Final time (min)
    
    dt=2;            % Time step size (min)
    if Stochastic==1
        dt=1;        % Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
    end
    Nt=(tf-t0)/dt;   % Number of time elements
    T=t0:dt:tf;      % Time vector
    
    VertProtrusions=0;         %Increases number of signalling neighbours in the vertical direction 
    HorzProtrusions=0;
    HorzProtStrength=1;
    
    %% Differentiation selection and parameters
    CrudeDiff           = 0;   %Cells will be marked...
    AvM                 = 0;   %Absolute vs moving mean threshold (0=Absolute thresh, 1=Moving mean thresh)
    wl_mins             = 100; %Window length in mins to take the moving mean from
    wl=wl_mins/dt;             %Window length converted to number of vector elements
    Replace             = 0;   %Replace differentiating cells with mean population protein
    DiffTime_hours=50;         %Time at which differention can start to occur (hours)
    DiffTime=DiffTime_hours*60/dt;
    S=0.01;                    % Rate of differentiation. Nominal value of 0.01
    
    %% Cell-movement/swapping parameters
    ImplementSwapping   = 0;   %Cells will randomly swap in the horizontal direction
    Pm=0.005; % Total probability of movement either left or right (max Pm value is 0.5)
    SwapThresh=1-Pm/2; % Treshold to be used in probability of swapping function within diffSolver.m
    Ts=5; % Time (in mins) between each swapping event (nominal Ts=5)       
    
   
    %% Passing model parameters to the solver
    modelParams.rows                  = rows;
    modelParams.cols                  = cols;
    modelParams.t0                    = t0;
    modelParams.tf_hours              = tf_hours;
    modelParams.Stochastic            = Stochastic;
    modelParams.CoupledCells          = CoupledCells;
    modelParams.Boundary              = Boundary;
    modelParams.TurnOffAutorepression = TurnOffAutorepression;
    modelParams.a_m                   = a_m;
    modelParams.a_p                   = a_p;
    modelParams.P_H0                  = P_H0;
    modelParams.TauH                  = TauH;
    modelParams.n_H                   = n_H;
    modelParams.u_m                   = u_m;
    modelParams.u_p                   = u_p;
    modelParams.P_ND0                 = P_ND0;
    modelParams.n_ND                  = n_ND;
    modelParams.TauND                 = TauND;
    modelParams.cells                 = cells;
    modelParams.tf                    = tf;
    modelParams.dt                    = dt;
    modelParams.Nt                    = Nt;
    modelParams.T                     = T;
    modelParams.VertProtrusions       = VertProtrusions;
    modelParams.HorzProtrusions       = HorzProtrusions;
    modelParams.HorzProtStrength      = HorzProtStrength;
    modelParams.percentCounter        = 0;
    

    %% Run model and plots
    degRate = [1.1; 1]; %This way around for plotting purposes
%     clear CoV_indCells CoV_avgIndCells;
    plotCols = 5;
    plotRows = ceil(NumberOfParamSets/plotCols);
    plotRows = 4;
    if plotRows > 5
        plotRows = 5;
    end

    if mod(paramSet, 50) == 0
        figNum = figNum + 1;
        subPlotNum = 1;
    end
    u_p = modelParams.u_p;
    for i = 1:2
        
        modelParams.u_p = u_p * degRate(i);
        
        [P, M, t_step, CT, ~] = runModel(modelParams);

        %Population CoV
        SD(i) = std(P(:, round(0.5*Nt):end), 0, 'all');
        AVG(i) = mean(P(:, round(0.5*Nt):end), 'all');
        Median(i) = mean(P(:, round(0.5*Nt):end), 'all');
        CoV(i) = SD(i)/AVG(i);

        %Individual cell CoV
        SD_indCells = std(P(:, round(0.5*Nt):end), 0, 2);
        AVG_indCells = mean(P(:, round(0.5*Nt):end), 2);
        CoV_indCells{i}(paramSet, :) = SD_indCells./AVG_indCells;
        CoV_avgIndCells(paramSet, i) = median(SD_indCells./AVG_indCells);

        %Plots
        figure(figNum)
        set(gcf,'renderer','painters') %For EPS file export
        subplot( plotRows, plotCols, subPlotNum)
        if i == 1
            plot(T/60, P, 'color', color2); hold on
        else
            plot(T/60, P, 'color', color1);
            title(sprintf('error=%.2f', CoV(2) - CoV(1)))
        end

        figure(figNum + 100)
        set(gcf,'renderer','painters') %For EPS file export
        subplot( plotRows, plotCols, subPlotNum)
        
        if i == 1
            histogram(P(:, 0.5*Nt:end), 'edgecolor', 'none', 'facecolor', color2(1:3), 'facealpha', 0.3,'Normalization', 'probability'); hold on
            yticks([])
        elseif i == 2
            histogram(P(:, 0.5*Nt:end), 'edgecolor', 'none', 'facecolor', color1(1:3), 'facealpha', 0.3, 'Normalization', 'probability');
            yticks([])
            title(sprintf('%.2f', COV(paramSet, 1)-COV(paramSet, 2)))
        end
    end

    subPlotNum = subPlotNum + 1;

end


%% Plot
CoV_avgIndCells = fliplr(CoV_avgIndCells); %This is the wrong way round to make a plot look nice in the previous section of code, so reverting it back here

figure(201)
subplot(131)
violinplot(CoV_avgIndCells,{'Normal protein deg', '1.1X degradation'}, 'QuartileStyle','shadow', 'ShowMean', true, 'ViolinColor', [color1(1:3); color2(1:3)]);
title('CoV')
ylabel('CoV')
axis square
drawnow
set(gca, 'fontsize', 15)

[pVal_CoV_indCells, ~] = kruskalwallis(CoV_avgIndCells); %Statistical test



