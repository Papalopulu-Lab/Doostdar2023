
function [err, coherence_mean, coherence_indCells, period_mean, period_indCells, CoV, abundance_mean]=objFunctionSummaryStats(X)


%% Define grid size
rows = X(11);         %Number of rows of cells in hexagonal grid
cols = X(12);         %Number of columns of cells in hexagonal grid

%% Define simulation time and step size
t0=0;            % Start time
tf_hours=100;    % Final time (hours)

%% Set up simulation inputs

%Important options!
Stochastic = X(13);              %1 = stochastic, 0 = deterministic  (deterministic currently has a bug) 
CoupledCells = X(14);            %Simulate with Notch-Delta coupling = 1, uncoupled = 0
Boundary = X(15);                %0 = 'hard' boundary, 1 = periodic boundary (set to 1)
autorepression = X(16);

if autorepression 
    TurnOffAutorepression = 0; %Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell = 1, with autorepression = 0
else
    TurnOffAutorepression = 1;
end
    
%% Model parameters


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



%% Other parameter set-up
cells=rows*cols; %Total number of cells
tf=tf_hours*60;  % Final time (min)

dt=2;            % Time step size (min)
if Stochastic==1
    dt=1;        % Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
end
Nt=(tf-t0)/dt;   % Number of time elements
T=t0:dt:tf;      % Time vector

VertProtrusions=0;         %Increases number of signalling neighbours in the vertical direction 
HorzProtrusions=0;
HorzProtStrength=1;

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
% modelParams.u_p                   = u_p;
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

%% Run model
Rd = [1 1.1]; % Rate of degredation
polyOrder=2;                  %Order of detrending polynomial in detrend.m function
frameTime=40;                 %Frame length in hours for detrending window
frameLength=frameTime*60./dt; %Conversion to window length in elements
frac=0.2;


coherence_mean     = zeros(1, 2);
coherence_indCells = zeros(cells, 2);
period_mean        = zeros(1, 2);
period_indCells    = zeros(cells, 2);
abundance_mean     = zeros(1,2);

CoV = zeros(1,2);

figure(100)
clf
for i = 1:2
    
    modelParams.u_p = Rd(i)*u_p;
    [P, ~, ~, ~, ~] = runModel(modelParams); %Run the model

    if i == 1
        subplot(211)
        plot(T/60, P, 'color', 'r'); hold on
    elseif i == 2
        subplot(212)
        plot(T/60, P, 'color', 'b');   
    end
    
    [~, ~, ~, ~, ~, avgIndCoh, ~, ~, ~, ~, ~, cohIndCells] = individualCoherence(T,P,polyOrder,frac,frameLength,Nt,dt);

    dominantPeriodWTSignificant = zeros(cells,1);
%     Fs = hours(1/60);
    for cell = 1:cells
        x = P(cell, frac*Nt:end);
        [~, ~, ~, domPerWTSig, ~] = significantCWT(x, dt, 1);
        dominantPeriodWTSignificant(cell) = hours(domPerWTSig);
    end
    % dominantPeriodWTSignificant
    %Summary model outputs
    coherence_mean(1,i)     = avgIndCoh;
    coherence_indCells(:,i) = cohIndCells; 
    period_mean(1, i)       = mean(dominantPeriodWTSignificant);
    period_indCells(:,i)    = dominantPeriodWTSignificant; 
    CoV(i)                  = std(P(:, round(0.5*Nt):end), 0, 'all')/mean(P(:, round(0.5*Nt):end), 'all');
    abundance_mean(1,i)     = mean(P(:, Nt/2:end),'all');
    
    
end

figure(100)
subplot(211)
legend(sprintf('CoV = %.2f', CoV(1)));
subplot(212)
legend(sprintf('CoV = %.2f', CoV(2)));

err = CoV(1) - CoV(2);
title(sprintf('Error = %.f', err))
drawnow

