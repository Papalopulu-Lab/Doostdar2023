%runOptimiser.m Runs a pattern search optimiser N-times starting from
%random initial parameter values, and returns a minimised error parameter
%set on each run. It optimises on a multicellular model of Her6 expression
%and minimises for maximal changes in population coefficient of variation
%when protein degradation rate is increased.

% REQUIRED: - Optimization toolbox
%           - Signal toolbox
%           - Statistics toolbox    
%           - Wavelet toolbox

clear;clc;
addpath('Functions')


%% How many times to run the optimiser
N = 6000;                             %Number of time to run the pattern search optimiser


%% Model setup parameters to pass to the objective function
rows        = 10;                     %Number of rows in the hexagonal grid of cells simulated
cols        = 6;                      %Number of columns in the hexagonal grid of cells simulated
stochastic  = true;                   %If stochastic == true, then simulations will be solved as chemical langevin equations using a Euler-Maruyama solver
coupled     = true;                   %If coupled == true, then cell Her6 dynamics will be coupled via lateral inhibition
autorepression = true;                %'true' means HES5 will repress its own expression.

if rows == 1 || cols == 1
    periodicBoundaries = false;       %Set non-periodic boundary conditions when simulation is 1D
else
    periodicBoundaries = true;        %Set periodic boundary conditions when simulation is 2D
end

modelSetup = [rows, cols, stochastic, coupled,  periodicBoundaries, autorepression]; %This is passed through the objective function optimisations parameters as it doesnt seem like there is a better way of doing it


%% Define function to be optimised
fun=@objFunction;                     %This function runs the model with the same parameter set twice, but the second run uses the 1.1X the degradation rate


%% |---------------------- Upper and lower bound parameter values -------------------------|     |-- Passing modelSetup to objective function here --|                                             
%   a_m   a_p    P_H0    n_H    Tau_H    P_ND0    n_ND    TauND   half-life m    half-life p            
LB=[0.1,   1,     100,    1,       0,      1,      1,        0,        2,            2,                           modelSetup]; %Lower parameter limits for optimisation
UB=[ 40,  40,   25000,    4,      20,    25000,    4,       60,       20,           20,                           modelSetup]; %Upper parameter limits for optimisation


%% Pattern search initilisation
XMIN = zeros(N, 10);                  %To store the value of all model parameters that minimises the optimiser error in each run of the optimiser
FVAL = zeros(N, 1);                   %To store all final error values that correspond to each minimised parameter set 
X_INIT=(UB-LB).*rand(N,numel(LB))+LB; %Generate random starting points in parameter space for each run of the pattern search method


%% Initialise arrays that will store model summary statistics for each identified parameter set
COHERENCE_MEAN = zeros(N, 2);         %To store the population mean coherence of temporal oscillations for each parameter set that the optimiser identifies
COHERENCE_IND_CELLS = cell(N,1);      %To store the coherence of temporal oscillations in individual cells
PERIOD_MEAN = zeros(N, 2);            %To store mean temporal period for each parameter set that the optimiser identifies
PERIOD_IND_CELLS = cell(N,1);         %To store the temporal period of individual cells
ABUNDANCE_MEAN = zeros(N, 2);         %To store mean population expression for each identified parameter set
COV = zeros(N, 2);                    %To store the population coefficient of variation for each parameter set


%% Other setup
dt = 1;                               %Time-step size in minutes


%% Run the pattern search method N times
for n = 1:N
    clc
    fprintf(sprintf('Optimiser run: %.f/%.f', n, N))
    
    X_init = X_INIT (n, :);           %Generate random starting point for patternsearch method
    
    PSoptions = optimoptions(@patternsearch, 'PlotFcn', {@psplotbestf,@psplotbestx}, 'UseParallel', true, 'UseCompletePoll', true, 'MaxFunctionEvaluations',500); %Set options for the pattern search optimiser
    [Xmin,fval,exitflag,output]=patternsearch(fun,X_init,[],[],[],[],LB,UB,PSoptions); %Run pattern search optimiser
    
    XMIN(n,:)=Xmin(1:10);
    FVAL(n)=fval;

    Fs = hours(dt/60); %Placing here because of a wierd bug with the wavelet transform function not liking 'duration' being converted to 'double' - this prevents that error, but I have no idea why!
    [err, coherence_mean, coherence_indCells, period_mean, period_indCells, CoV, abundance_mean] = objFunctionSummaryStats([Xmin, modelSetup]); %Run outside of patternsearch with the final parameter set to return useful values for filtering the data further.
    
    %Store model summary statistics
    COHERENCE_MEAN(n, :) = coherence_mean;
    COHERENCE_IND_CELLS{n} = coherence_indCells;
    PERIOD_MEAN(n,:) = period_mean;
    PERIOD_IND_CELLS{n} = period_indCells;
    ABUNDANCE_MEAN(n, :) = abundance_mean;
    COV(n, :) = CoV;
    
    %Select values that have error below a certain value to plot
    cond2 = (max(COHERENCE_MEAN(1:n, :), [], 2) < 0.5);

    elemKeep=find(cond2 == 1);
    XMINplot=XMIN(elemKeep, :);
    FVALplot=FVAL(elemKeep);
    
    fprintf(sprintf('\n Percentage that are below error threshold = %.f', 100*length(elemKeep)/N))
    
    axisLabels={'a_m', 'a_p', 'P_{H0}', 'n_H', 'Tau_H', 'P_{ND0}', 'n_{ND}', 'Tau{ND}', 'mRNA T_{1/2} ', 'protein T_{1/2}'};
    figure(9)
    clf
    for m=1:10
        subplot(2,5,m)
        histogram(XMINplot(:,m))
        xlim([LB(m) UB(m)])
        axis square
        xlabel(axisLabels(m))
    end

    %Scatter plot
    figure(10)
    param1 = 10;
    param2 = 3;
    param3 = 6;
    scatter3(XMINplot(:, param1), XMINplot(:, param2), XMINplot(:, param3), 40, FVALplot, 'filled')
    xlabel('Protein half-life (min)')
    ylabel('P_{H0}')
    zlabel('P_{ND0}')
    colorbar
    drawnow
    
end


%% Save optimiser output
STR = sprintf("Data/Autorepression=%.0f Coupled=%.0f Rows=%.0f  Cols=%.0f  Stochastic=%.0f  Number of runs = %.0f.mat", autorepression, coupled, rows, cols, stochastic, N);
save(STR, 'XMIN', 'FVAL', 'COHERENCE_MEAN', 'COHERENCE_IND_CELLS', 'PERIOD_MEAN', 'PERIOD_IND_CELLS', 'ABUNDANCE_MEAN', 'COV', 'LB', 'UB', 'modelSetup');


