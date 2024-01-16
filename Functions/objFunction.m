
function [err] = objFunction(X)


%% Define grid size
rows = X(11);                  %Number of rows of cells in hexagonal grid
cols = X(12);                  %Number of columns of cells in hexagonal grid

%% Define simulation time and step size
t0=0;                          % Start time
tf_hours=100;                  % Final time (hours)

%% Set up simulation inputs

Stochastic     = X(13);        %1 = stochastic, 0 = deterministic  (deterministic currently has a bug) 
CoupledCells   = X(14);        %Simulate with Notch-Delta coupling = 1, uncoupled = 0
Boundary       = X(15);        %0 = 'hard' boundary, 1 = periodic boundary (set to 1)
autorepression = X(16);

if autorepression 
    TurnOffAutorepression = 0; %Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell = 1, with autorepression = 0
else
    TurnOffAutorepression = 1;
end
%% Model parameters


a_m   = X(1);                  %mRNA production rate
a_p   = X(2);                  %Protein production rate
P_H0  = X(3);                  %Autorepression repression threshold
n_H   = X(4);                  %Autorepression Hill coefficient
TauH  = X(5);                  %Autorepression time delay

P_ND0 = X(6);                  %Lateral inhibition repression threshold
n_ND  = X(7);                  %Lateral inhibition Hill coefficient
TauND = X(8);                  %Lateral inhibition time delay

u_m   = log(2)/X(9);           %mRNA degredation rate
u_p   = log(2)/X(10);          %Protein degredation rate



%% Other parameter set-up
cells=rows*cols;               %Total number of cells
tf=tf_hours*60;                %Final time (min)

dt=2;                          %Time step size (min)
if Stochastic==1
    dt=1;                      %Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
end
Nt=(tf-t0)/dt;                 %Number of time elements
T=t0:dt:tf;                    %Time vector

VertProtrusions=0;             %Increases number of signalling neighbours in the vertical direction 
HorzProtrusions=0;             %Increases number of signalling neighbours in the horizontal direction 
HorzProtStrength=1;            %Horizontal protrusion coupling strength

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
Rd = [1 1.1];       % Rate of degredation
CoV = zeros(1,2);   % Coefficient of variation

for i = 1:2
    
    modelParams.u_p = Rd(i)*u_p;
    [P, ~, ~, ~, ~] = runModel(modelParams); %Run the model
    CoV(i)          = std(P(:, round(0.5*Nt):end), 0, 'all')/mean(P(:, round(0.5*Nt):end), 'all');
  
end
err = CoV(1) - CoV(2);


