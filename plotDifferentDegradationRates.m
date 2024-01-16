% plotDifferentDegradationRates.m Uses one parameter set from the optmiser
% parameter sets in the Data folder and simulates the same parameter set
% over several different protein degradation rates, and then plots these
% outputs as histograms of Her6 expression distributions.

clear;clc;
addpath('Functions')
addpath('Data')


%% Choose the colour of all graphs
GraphAppearance = 0;          %0 = Normal, 1 = Dark grey

if GraphAppearance==1   
    BackgroundColour=38/255*[1 1 1]; TextColour=[1 1 1];
    get(0,'Factory');                           set(0,'defaultfigurecolor',BackgroundColour)
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor',BackgroundColour)
    set(0,'defaultAxesXColor',TextColour);      set(0,'defaultAxesYColor',TextColour)
    set(0,'defaultLegendTextColor',TextColour); set(0,'defaultTextColor',TextColour)  
    set(0,'defaultAxesZColor',TextColour)
elseif GraphAppearance==0
    BackgroundColour = [1 1 1];
    get(0,'Factory');                           set(0,'defaultfigurecolor',[0.94 0.94 0.94])
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor','w')
    set(0,'defaultAxesXColor','k');             set(0,'defaultAxesYColor','k')
    set(0,'defaultLegendTextColor','k');        set(0,'defaultTextColor','k')
end


%% Define grid size
rows=10;                     %Number of rows of cells in hexagonal grid
cols=6;                      %Number of columns of cells in hexagonal grid


%% Define simulation time and step size
t0=0;                        %Start time
tf_hours=100;                %Final time (hours)


%% Set up simulation inputs
%Important options!
Stochastic = 1;              %1 = stochastic, 0 = deterministic  (deterministic currently has a bug) 
CoupledCells = 1;            %Simulate with Notch-Delta coupling = 1, uncoupled = 0
Boundary = 1;                %0 = 'hard' boundary, 1 = periodic boundary (set to 1)
TurnOffAutorepression = 0;   %Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell = 1, with autorepression = 0


%% Other parameter set-up
cells=rows*cols;            %Total number of cells
tf=tf_hours*60;             % Final time (min)

dt=2;                       % Time step size (min)
if Stochastic==1
    dt=1;                   % Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
end
Nt=(tf-t0)/dt;              % Number of time elements
T=t0:dt:tf;                 % Time vector

VertProtrusions=0;          %Increases number of signalling neighbours in the vertical direction 
HorzProtrusions=0;
HorzProtStrength=1;


%% Filter data
load("Model 2 Coupled - Rows=10  Cols=6  Stochastic=1  Number of runs = 6000.mat"); % Model 2: Coupled parameter set - remember to set the 'CoupledCells' variable above to 1 if using these parameter sets

fval_cutoff      = -0.25;   %Filter for values less than this (nominal = -0.25)
coherence_cutoff = 0.5;     %Filter for values less than this (nominal = 0.5)
abundance_cutoff = 2000;    %Filter for values greater than this (nominal = 3000)
elemKeep=find(FVAL < fval_cutoff   &   max(COHERENCE_MEAN, [], 2) < coherence_cutoff    &   ABUNDANCE_MEAN(:,1) > abundance_cutoff == true);

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
figNum = 1;
subPlotNum = 1;
modelParams={};

X = XMIN(63,:);

%% Loop over different degradation rates
II=4; 
JJ=1; 
ii_arr=linspace(1,1.1,II);%Degredation rate
jj_arr=[linspace(1000,1300,JJ-1) 100000];%Coupling strength
P_cell = cell(II, JJ);

if JJ==1
    jj_arr=X(6);
end
steadyState=zeros(cells,II,JJ);
count=0;

for ii = 1:II %Degradation rate
        
    count=count+1;
    clc
    fprintf(sprintf('Progress: %.f/%.f \n',count,II*JJ))
    a_m   = X(1);         %mRNA production rate
    a_p   = X(2);         %Protein production rate
    P_H0  = X(3);         %Autorepression repression threshold
    n_H   = X(4);         %Autorepression Hill coefficient
    TauH  = X(5);         %Autorepression time delay
    
    P_ND0 = X(6);         %Lateral inhibition repression threshold
    n_ND  = X(7);         %Lateral inhibition Hil coefficient
    TauND = X(8);  
    
    u_m   = log(2)/X(9);  %mRNA degredation rate
    u_p   = ii_arr(ii)*log(2)/X(10); %Protein degredation rate

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
    modelParams.percentCounter        = 1;
        

    %% Run model
    [P, M, t_step, CT, ~] = runModel(modelParams);
    P_cell{ii} = P;
end


%% Plots
figure(1211)
set(gcf,'renderer','painters') %For EPS file export
clf
load cmap


plotCols = 10;
bigSubplot = 1:II*plotCols;
bigSubplot(10:10:end) = [];

ha=subplot(II,10, [bigSubplot]);

col1 = [77 66 66]./255;
col2 = [229, 98, 129]./255;
col3 = [252 238 242]./255;

colour = colourMapInterp(col1, col2, col3, 100);

figure(1222)
clf
for ii = 1:II
    d = P_cell{ii}(:, 0.5*Nt:end);
    data(ii, :) = d(:); 
    
    subplot(II,1, ii)
    histogram(d, 50, 'FaceColor', colour(100*(ii-1)/II +1,:),'normalization' , 'probability')
    xlim([0 26000])
    ylim([0 0.2])
    set(gca, 'fontsize', 20)
end
    
figure(1211)

    edges = linspace(min(data,[],'all'), max(data,[],'all'), 100); % bin edges
    counts = histc(data, edges, 2); % specify dim 2 to act column-wise

    fillSize=15; %Size of gap between 

    fill=zeros(size(counts,1)*fillSize, size(counts,2));
    customCmap=repmat(BackgroundColour,size(counts,1)*fillSize,1);
    colourMap=cmap;
    colourMap=randColourMapBright(100,[0 0 0]);
    countElem=0;

    for i=1:fillSize:size(fill,1)
        countElem=countElem+1;
        fill(i,:)=counts(countElem,:);
    end
    fill(end - fillSize + 2:end,:) = [];
    counts = fill;
    
    hb = bar3(edges, counts.',1); % note the transpose to get the colors right
    hbh = get(hb(3),'parent');

    xlabel('Degradation factor')
    ylabel('Her6 abundance')
    ylim([edges(1)-(edges(end)-edges(end-1)) edges(end)+(edges(end)-edges(end-1))])
    zlabel('Count');

    colormap([BackgroundColour; flipud(cmap)])
    view(320,40)
    set(gca,'FontSize',11, 'ZColor', 'k')
    zlim([0 max(counts,[],'all')])
    set(hb,'EdgeColor','none');

    numXTicks = 5;
    xticks(linspace(1,size(fill,1),numXTicks))
    xticklabels(num2str(linspace(ii_arr(1), ii_arr(end), numXTicks)','%.2f'))

    colour = flipud(brewermap(100,'Greens'));
    col1 = 0.2*[1 1 1];
    col2 = 1/255*[156, 213, 153];

    colour = gray(100);

    colour=colour(1:80,:);
    colormap(colour); 

for i = 1:numel(hb)
  index = logical(kron(counts(i,:) == 0, ones(6, 1)));
  zData = get(hb(i), 'ZData');
  zData(index, :) = nan;
  set(hb(i), 'ZData', zData);
end

view([269.99 78])
set(gca,'FontSize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

for ii = 1:II
    minAbundance(ii) = min(P_cell{ii,1}(:, Nt/2:end), [], 'all');
    maxAbundance(ii) = max(P_cell{ii,1}(:, Nt/2:end), [], 'all');
end
minAbundance = min(minAbundance);
maxAbundance = max(maxAbundance);

for ii = 1:II
figure(1211)
    
    fig=gcf;
    if GraphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
    colour_index=reshape(flipud(vecTOmat(P_cell{ii}(:,t_step),cols)),[1,cols*rows]);

    map=magma(1000);

    ha = subplot(II, plotCols, (II+1-ii)*plotCols);
    hexagons = patch(X_vertices,Y_vertices,colour_index,'edgecolor','none');
    colormap(ha,map); 
    axis equal;
    set(gca,'Visible','off')
    caxis([minAbundance 0.4*maxAbundance])
    colour_index=reshape(flipud(vecTOmat(P_cell{ii}(:,end),cols)),[1,cols*rows])';
    set(hexagons, 'FaceVertexCData',colour_index);

    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'FontSize',10)

    drawnow


end