% 
function [P, M, t_step, CT, HILL] = runModel(modelParams)

    rows                  = modelParams.rows;
    cols                  = modelParams.cols;
    t0                    = modelParams.t0;
    tf_hours              = modelParams.tf_hours;
    Stochastic            = modelParams.Stochastic;
    CoupledCells          = modelParams.CoupledCells;
    Boundary              = modelParams.Boundary;             
    TurnOffAutorepression = modelParams.TurnOffAutorepression;
    a_m                   = modelParams.a_m;
    a_p                   = modelParams.a_p;
    P_H0                  = modelParams.P_H0;
    TauH                  = modelParams.TauH;
    n_H                   = modelParams.n_H;
    u_m                   = modelParams.u_m;
    u_p                   = modelParams.u_p;
    P_ND0                 = modelParams.P_ND0;
    n_ND                  = modelParams.n_ND;
    TauND                 = modelParams.TauND;
    percentCounter        = modelParams.percentCounter;
    
    cells                 = modelParams.cells;
    tf                    = modelParams.tf;
    dt                    = modelParams.dt;
    Nt                    = modelParams.Nt;
    T                     = modelParams.T;
    VertProtrusions       = modelParams.VertProtrusions;
    HorzProtrusions       = modelParams.HorzProtrusions;
    HorzProtStrength      = modelParams.HorzProtStrength;





    %% Other parameter set-up

%     cells=rows*cols; %Total number of cells
%     tf=tf_hours*60;  % Final time (min)

%     dt=2;            % Time step size (min)
%     if Stochastic==1
%         dt=1;        % Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
%     end
%     Nt=(tf-t0)/dt;   % Number of time elements
%     T=t0:dt:tf;      % Time vector

%     VertProtrusions=0;         %Increases number of signalling neighbours in the vertical direction 
%     HorzProtrusions=0;
%     HorzProtStrength=1;

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


    %% Produce neighbour matrix

    if cells==1
        Boundary=0; % Other boundaries don't make sense for 1 cell!
    end

    [NM,NumNeigh]=neighbours(rows, cols, Boundary, VertProtrusions, HorzProtrusions, HorzProtStrength); 
    NM=sparse(NM); % Saves a lot on computational cost for large grid sizes!

    eps=(1./NumNeigh); % Make this as an output of neighbours.m!
    if cells==1
        eps=0;
    end


    %Experimental changing upper and lower bounds of intercellular Hill function
    Up=1;    % Upper limit of Hill function (nominally = 1)   
    Low=0.0; % Lower limit (nominally 0)

    %Experimental Cis-inhibition stuff - set min_Up to 1 for no cis-inhibition
    P_0Up=2e4;
    n_Up=2;
    min_Up=1; %Set this to one to make cis-inhibition not dependent upon Hes levels


    %__________________________________________________________________________

    gamma=1; % Maximum intercellular Hill function value

    TauH_step  = round(TauH/dt);  % Conversion to simulation time steps
    TauND_step = round(TauND/dt); % Conversion to simulation time steps

    if CoupledCells==0
        eps=0; %Set eps=0 to decouple the cells 
    end

    if CrudeDiff==0
        DiffTime=Nt*2;
    end



    %% DDEs (using anonymous functions (@ symbol) rather than standard functions)
    dm=@(m,p,p_delay,Psum_delay,gamma,a,b) a_m*1./(1+(p_delay./P_H0).^n_H).*gamma.*(  a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)  ) - u_m.*m; % Describe mRNA change in time


    if TurnOffAutorepression==1
        dm=@(m,p,p_delay,Psum_delay,gamma,a,b) a_m.*gamma.* (a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)) - u_m.*m; % No self repression version for simple Notch Delta modelling
    end

    dp=@(m,p,p_delay,Psum_delay,gamma,a,b) a_p.*m - u_p.*p; % Describes protein change in time

    %% SDDEs
    dm1=dm;

    dm2=@(m,p,p_delay,Psum_delay,gamma,a,b) sqrt(a_m.*gamma./(1+(p_delay./P_H0).^n_H).*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);

    if TurnOffAutorepression==1
       dm2=@(m,p,p_delay,Psum_delay,gamma,a,b) sqrt(a_m.*gamma.*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);
    end

    dp1=dp;
    dp2=@(m,p,p_delay,Psum_delay,gamma,a,b) sqrt(a_p.*m + u_p.*p);


    %% Initialise vectors and initial values (random)

    rng(3); % Same seed for random number generator (prevents initial heterogeneities from changing on every simulation, so comment out this line if you want different initial conditions each time)
    P=[rndrng(cells,1,5000,20000),zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation
    M=[rndrng(cells,1,0,20),zeros(cells,Nt-1),zeros(cells,1)];              %Vector that will store mRNA values from the simulation
    rng('default')
%     P=[linspace(5000, 20000, cells)',zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation
%     M=[linspace(0, 20, cells)',zeros(cells,Nt-1),zeros(cells,1)];              %Vector that will store mRNA values from the simulation

    %% Main loop (solves differential equations)
%     percentCounter=1; %Print progress of the differential solver to the command line (1=yes, 0=no)

    solverParams.P                 = P;
    solverParams.M                 = M;
    solverParams.Nt                = Nt;
    solverParams.TauND_step        = TauND_step; 
    solverParams.TauH_step         = TauH_step;
    solverParams.NM                = NM;
    solverParams.gamma             = gamma;
    solverParams.dm                = dm;
    solverParams.dp                = dp;
    solverParams.dm1               = dm1;
    solverParams.dm2               = dm2; 
    solverParams.dp1               = dp1; 
    solverParams.dp2               = dp2;
    solverParams.dt                = dt;
    solverParams.Stochastic        = Stochastic;
    solverParams.rows              = rows;
    solverParams.cols              = cols;
    solverParams.DiffTime          = DiffTime;
    solverParams.S                 = S;
    solverParams.ImplementSwapping = ImplementSwapping;
    solverParams.SwapThresh        = SwapThresh;
    solverParams.Ts                = Ts;
    solverParams.percentCounter    = percentCounter;
    solverParams.AvM               = AvM;
    solverParams.Replace           = Replace;
    solverParams.wl                = wl;
    solverParams.Up                = Up;
    solverParams.Low               = Low;
    solverParams.P_ND0             = P_ND0;
    solverParams.n_ND              = n_ND;
    solverParams.P_0Up             = P_0Up;
    solverParams.n_Up              = n_Up;
    solverParams.min_Up            = min_Up;
    

    [P, M, t_step, CT,HILL]=diffSolver(solverParams);

%     longExposure=DiffYNflash;
%     DiffElems=find(DiffYNflash==1);
%     ExposureTime=50;
%     for E=1:ExposureTime
%         Elems=DiffElems+cells*E;
%         Elems(Elems>numel(longExposure))=[];
%         longExposure(Elems)=1;
%     end


end




