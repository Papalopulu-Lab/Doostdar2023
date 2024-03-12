%runSingleParameterSet.m Uses a single parameter set from the
%optimiser-identified parameter sets in the Data folder, and plots the
%outputs of a single run of the model.

clear;clc;
addpath('Functions')
addpath('Data')

%% Choose the colour of all graphs
GraphAppearance = 0;                            %0 = Normal, 1 = Dark grey

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
rows=10;                                        %Number of rows of cells in hexagonal grid
cols=6;                                         %Number of columns of cells in hexagonal grid


%% Define simulation time and step size
t0=0;                                           %Start time
tf_hours=100;                                   %Final time (hours)


%% Set up simulation inputs
%Important options!
Stochastic = 1;                                 %1 = stochastic, 0 = deterministic  (deterministic currently has a bug) 
CoupledCells = 1;                               %Simulate with Notch-Delta coupling = 1, uncoupled = 0
Boundary = 1;                                   %0 = 'hard' boundary, 1 = periodic boundary (set to 1)
TurnOffAutorepression = 0;                      %Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell = 1, with autorepression = 0


%% Filter data
load("Model 2 Coupled - Rows=10  Cols=6  Stochastic=1  Number of runs = 6000.mat"); % Model 2: Coupled parameter set - remember to set the 'CoupledCells' variable above to 1 if using these parameter sets

fval_cutoff      = -0.25;                       %Filter for values less than this (nominal = -0.25)
coherence_cutoff = 1;                           %Filter for values less than this (nominal = 0.5)
abundance_cutoff = 2000;                        %Filter for values greater than this (nominal = 3000)
elemKeep=find(FVAL < fval_cutoff   &   max(COHERENCE_MEAN, [], 2) < coherence_cutoff    &   ABUNDANCE_MEAN(:,1) > abundance_cutoff == true);

XMIN=XMIN(elemKeep, :);                         %The filtered set of minimised parameter sets 
FVAL=FVAL(elemKeep);                            %The filtered set of optimiser error values 
COHERENCE_MEAN =  COHERENCE_MEAN(elemKeep, :);  %Mean coherence of each filtered parameter set
PERIOD_MEAN = PERIOD_MEAN(elemKeep, :);         %Mean period of each filtered parameter set
COV = COV(elemKeep, :);                         %Coefficient of variation for each filtered parameter set

for elem = 1 : numel(elemKeep)
    COHERENCE_IND_CELLS_Filtered{elem, 1} = COHERENCE_IND_CELLS{elemKeep(elem)}; %Coherence of individual cells within each parameter set simulation
    PERIOD_IND_CELLS_Filtered{elem, 1} = PERIOD_IND_CELLS{elemKeep(elem)};       %Period of individual cells within each parameter set simulation
end


%% Order parameter sets by a certain output value
paramSet = 47;                                  %Choose which parameter set to simulate and plot

[FVAL,sortIdx] = sort(FVAL);                    %Reorder by largest FVAL 

XMIN = XMIN(sortIdx,:);                         %Reorder by largest FVAL
FVAL = FVAL(sortIdx);                           %Reorder by largest FVAL
COHERENCE_MEAN = COHERENCE_MEAN(sortIdx, :);    %Reorder by largest FVAL
PERIOD_MEAN = PERIOD_MEAN(sortIdx, :);          %Reorder by largest FVAL
COV = COV(sortIdx, :);                          %Reorder by largest FVAL

expectedCOV = COV(paramSet, :);              

for idx = 1 : numel(sortIdx)
    COHERENCE_IND_CELLS_Filtered{idx, 1} = COHERENCE_IND_CELLS{sortIdx(idx)};
    PERIOD_IND_CELLS_Filtered{idx, 1} = PERIOD_IND_CELLS{sortIdx(idx)};
end

NumberOfParamSets = length(FVAL);
figNum = 1;
subPlotNum = 1;
modelParams={};

X = XMIN(paramSet, :);

a_m   = X(1);                                   %mRNA production rate
a_p   = X(2);                                   %Protein production rate
P_H0  = X(3);                                   %Autorepression repression threshold
n_H   = X(4);                                   %Autorepression Hill coefficient
TauH  = X(5);                                   %Autorepression time delay

P_ND0 = X(6);                                   %Lateral inhibition repression threshold
n_ND  = X(7);                                   %Lateral inhibition Hill coefficient
TauND = X(8);                                   %Lateral inhibition time delay

u_m   = log(2)/X(9);                            %mRNA degredation rate
u_p   = 1*log(2)/X(10);                         %Protein degredation rate


%% Graphing/visualising output options
%Various frequency analysis options
TemporalFourier     = 1;                        %Gives average Fourier period of the cell population
TemporalWavelet     = 1;                        %Preliminary implementation of wavelet to examine whether Her6 switches between periodic and noisy

%Animation
AnimateGrid         = 0;                        %Set as 1 to animate the Her6 expression
AnimationSpeed      = .1;                       %Multiplier for speed of animation
ShowRandomCells     = 1;                        %Time traces of individual cells
ShowLastFrame       = 1;                        %Shows the last time point in the hexagonal lattice arrangement of cells

%Saving animation options
MakeGIF   = 0;                                  %1=yes, 0=no. Any animations that run will be made into GIFs.
filename1 = 'grid animation';                   %Specify the output file name


%% Other parameter set-up
cells=rows*cols;                                %Total number of cells
tf=tf_hours*60;                                 %Final time (min)

dt=2;                                           %Time step size (min)
if Stochastic==1
    dt=1;                                       %Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
end
Nt=(tf-t0)/dt;                                  %Number of time elements
T=t0:dt:tf;                                     %Time vector

VertProtrusions=0;                              %Increases number of signalling neighbours in the vertical direction 
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
[P, M, t_step, CT, ~,] = runModel(modelParams);

CoV = std(P(:, round(0.5*Nt):end), 0, 'all')/mean(P(:, round(0.5*Nt):end), 'all');
fprintf(sprintf('Coefficient of variation for the population = %.2f', CoV))


%% Make table of parameters used and print to command line
Parameter = ["a_m"; "a_p"; "P_H0"; "n_H"; "TauH"; "P_ND0"; "n_ND"; "TauND"; "u_m"; "u_p"];
% Age = [38;43;38;40;49];
Value = X';

Parameters = table(Parameter,Value)


%==========================================================================
%%                        Animation of Grid
%==========================================================================

if AnimateGrid==1
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
    
    colour_index1=reshape(flipud(vecTOmat(P(:,1),cols)),[1,cols*rows]);
    
    figure(1)
    clf;
    fig=gcf;
    fig.InvertHardcopy = 'on';

    map=viridis(1000);
    colormap(map);

    hexagons1 = patch(X_vertices,Y_vertices,colour_index1,'edgecolor','none');
    set(gca,'xtick',[],'ytick',[]); 

    axis equal; title('Hes expression')
    set(gca,'Visible','off')
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'fontsize',15)
    caxis([min(min(P(:,round(0.2*Nt):end))) max(max(P(:,round(0.2*Nt):end)))])
    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    
    startT=Nt*0.5;
    T_STEP=startT:round(AnimationSpeed*Nt/400):Nt;
    T_STEP=1:round(AnimationSpeed*Nt/400):Nt;
    ti=length(T_STEP);
    im=struct([]);
    
    for idx=1:ti
        
        t_step=T_STEP(idx);
        
        colour_index1=reshape(flipud(vecTOmat(P(:,t_step),cols)),[1,cols*rows])';
        set(hexagons1, 'FaceVertexCData',colour_index1); 
        title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
        drawnow;  

        if MakeGIF==1
            F1=getframe(gcf); %gca makes frame without labels, gcf makes frame including labels
            im{idx}=frame2im(F1);
            F_all{idx}=F1;
        end
    end
    
    if MakeGIF==1
        IDX=idx;
        for idx = 1:IDX
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,filename1,'gif','LoopCount',Inf,'DelayTime',1/30);
            else
                imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',1/30);
            end
        end
        
        video = VideoWriter(filename1); %create the video object
        video.FrameRate = 30;
        open(video); %open the file for writing
        for idx = 1:IDX
        %     I = imread(im{idx}); %read the next image
            writeVideo(video,F_all{idx}); %write the image to file
        end
        close(video); %close the file
    end
end


%==========================================================================
%%              Final time point of hexagonal grid plot
%==========================================================================
if ShowLastFrame==1
    figure(103)
    clf;
    
    fig=gcf;
    if GraphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
    colour_index=flipud(TauH);
    map=viridis(1000);
    hexagons = patch(X_vertices,Y_vertices,colour_index,'edgecolor','none');
    colormap(map); 
    colorbar

    set(gca,'Visible','off')  
    colour_index=reshape(flipud(vecTOmat(P(:,end),cols)),[1,cols*rows])';
    set(hexagons, 'FaceVertexCData',colour_index);
    
    title('Her6 abundance')
    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'FontSize',20)
    axis equal

    drawnow
end


%==========================================================================
%%                Plots of Hes levels in random cells 
%==========================================================================
if ShowRandomCells==1
    
    % Figure preamble
    figure(2) 
    clf
    set(gcf,'renderer','Painters')
    fig=gcf;
    if GraphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    
    
    %How many cells to plot
    desired_cells=cells;
    desired_fraction=desired_cells/cells;
    rand_ind=ceil(cells*rand(1,ceil(cells*desired_fraction)));
    
    if cells==2
        rand_ind=[1 2];
    end
    tStart=Nt-20*60/dt;
    num=length(rand_ind);
    subplot(1,5,[1 2 3 4])
    h=plot(T(tStart:end)/60,P(rand_ind,tStart:end),'linewidth',1); hold on
    set(h, {'color'}, num2cell(randColourMap(num),2));
    xlabel('Time (hours)')
    ylabel(sprintf('Hes protein count'))
    xlim([tStart*dt/60 tf_hours])
    set(gca,'FontSize',15)
    title('Subset of cell Her6 time traces')
    ylim([0 1.01*max(max(P(:,tStart:end)))])
    
    subplot(1,5,5)
    histogram(P(:,Nt/2:Nt),'normalization','probability','edgecolor','none','facecolor',[0.4 0.4 0.4])
    xlim([0 1.01*max(max(P(:,Nt/2:Nt)))])
    ylabel('Density')
    title(sprintf('Median=%.f',median(P(:,Nt/2:Nt),'all')))
    set(gca,'view',[90 -90]) 
end


%==========================================================================
%%                          Fourier transform
%==========================================================================
    
    %Changes fraction of signal to use (deterministic uses earlier part of
    %the signal as it is a damped oscillator)
    if Stochastic==1
        frac=0.2;
    else   
        frac=0.01;
    end
    Osc=find(max(P(:,round(0.7*Nt):Nt),[],2)>100);
    Y_raw=P(Osc,:);
     
    %Detrending
    polyOrder=2;                  %Order of detrending polynomial in detrend.m function
    frameTime=40;                 %Frame length in hours for detrending window
    frameLength=frameTime*60./dt; %Conversion to window length in elements  
%     frameLength=75;
    [Ydetrend,t,Ysmooth,f,P1,coherence,f_C1,P_C1,avgFourier,ind_per,I]=tempFourier(T,Y_raw,polyOrder,frac,frameLength,Nt,dt);
    
    
if TemporalFourier==1   
      
    fprintf('Coherence of oscillators in this model = %.2f\n',coherence)
    fprintf('Oscillators with periods below 7 hours = %.2f%%\n', 100*sum(ind_per<7)/numel(ind_per))
   
% Plots
    color1 = [0.55 0.8 0.65];
    color2 = [0.9 0.9 0.9];
    
    figure(14)
    clf
    fig=gcf;
    fig.InvertHardcopy = 'off';
    txtSize=10;
    sampleCell=size(Y_raw,1);
    
    subplot(321)
    plot(T(frac*Nt:end)./60,Y_raw(sampleCell,frac*Nt:end),'linewidth',1,'color', color1); hold on;
    plot(t./60,Ysmooth(sampleCell,:),'--', 'color','k');
    title('Signal with detrend line')
    xlim([dt*Nt*frac/60 dt*Nt/60])
    xlabel('Time (hours)')
    set(gca,'Fontsize',txtSize)
    
    subplot(322)
    plot([t(1)/60 t(end)/60], [0 0],'k'); hold on
    plot(t./60,Ydetrend(sampleCell,:),'color', color1,'linewidth',1); 
    xlim([t(1)./60 t(end)./60])
    title('Detrended signal')
    xlabel('Time (hours)')
    set(gca,'Fontsize',txtSize)
    
    subplot(3,2,3)
    h0=plot((60*f), P1, 'color', [0 0 0 0.1]);
    grid on
    xlim([0 0.5])
    xlabel('Frequency (1/h)')
    ylim([0 1.2*max(max(P1))])
    title('Individual cell Fourier transforms')
    set(gca,'Fontsize',txtSize)
    
    subplot(324)
    histo=histogram(ind_per,'normalization','probability');
    set(histo,'edgecolor','k');
    set(histo,'facecolor', color2);
    xlabel('Dominant period (h)')
    ylabel('Num of cells')
    title('Individual cell dominant periods')
    set(gca,'fontsize',txtSize)
    
    subplot(3,2,[5 6])
    h1=area(60*f,avgFourier); hold on;
    h1.FaceColor=[0.95 0.95 0.95];
    h2=area(60*f_C1, P_C1); hold on;
    h2.FaceColor = color1;
    h2.FaceAlpha = 0.3;
    plot(60*f,avgFourier,'.','linewidth',1,'color','k'); 
    plot(60*f(I),avgFourier(I),'.','markersize',20,'color', [0 0 0])
    title('All cells Fourier transform')
    xlabel('Frequency (1/h)')
    xlim([0 0.5])
    ylim([0 1.2*max(avgFourier)])
    title(sprintf('Sum of Fourier Transforms   |   Dominant period: %.2f hours   |   Coherence: %.2f', 1./(60*f(I)), coherence))
    set(gca,'Fontsize',txtSize)
    drawnow
 
%__________________________________________________________________________
   
    numPlots=24;
    if numPlots>cells
        numPlots=cells;
    end
    
    rand_ind=randperm(size(Y_raw,1));
    rand_ind=rand_ind(1:numPlots);
    
    t_start_h=tf_hours-40;
    figure(88)
    clf
    fig=gcf;
    fig.InvertHardcopy = 'off';
    
    for n=1:numPlots
      
        subplot(4,6,n) 
         
        h12=plot(T(:)/60,Y_raw(rand_ind(n),:),'linewidth',1); hold on
        plot(t./60, Ysmooth(rand_ind(n),:),'k--');
        set(h12, {'color'}, num2cell(randColourMapBright(1,[0.2 0.5 0.2]),2));
        xlabel('Time (hours)')
        ylabel(sprintf('Hes protein'))
        set(gca,'FontSize',8)
        xlim([t_start_h, tf_hours])
        title(sprintf('T=%.1fh',ind_per(rand_ind(n))))
        
    end
    drawnow

    figure(89)
    clf
    fig=gcf;
    fig.InvertHardcopy = 'off';
    
    for n=1:numPlots
      
        subplot(4,6,n) 
        h12=plot(t./60,Ydetrend(rand_ind(n),:),'linewidth',1); hold on
        plot(t./60, Ysmooth(rand_ind(n),:)*0,'k--');
        set(h12, {'color'}, num2cell(randColourMapBright(1,[0.2 0.5 0.2]),2));
        xlabel('Time (hours)')
        ylabel(sprintf('Hes protein'))
        set(gca,'FontSize',8)
        xlim([t_start_h, tf_hours])
        title(sprintf('T=%.1fh',ind_per(rand_ind(n))))
        
    end
    drawnow

    polyOrder=2;                  %Order of detrending polynomial in detrend.m function
    frameTime=40;                 %Frame length in hours for detrending window
    frameLength=frameTime*60./dt; %Conversion to window length in elements
    frac=0.2;
    [~, ~, ~, ~, ~, avgIndividualCoherence, ~, ~, ind_per, FisherPer, Occurrence, coherenceIndividualCells] = individualCoherence(T,P,polyOrder,frac,frameLength,Nt,dt);
    
    figure(15)
    clf
    % histogram(coherenceIndividualCells)
    meanExpression = mean(P(:, 0.2*Nt:end), 2);
    scatter(coherenceIndividualCells, meanExpression)
    title('Coherence versus expression level')
    ylabel('Abundance')
    xlabel('Coherence')
end

%==========================================================================
%%                          Wavelet Transform
%==========================================================================

if TemporalWavelet==1
    
x = P(5, 0.2*Nt:end);
t_WT = T(0.2*Nt:end);
randIdx = randperm(length(x));
xRand = x(randIdx);

numRepeats = 3;
[WT, dominantPeriodWT, significantWT, dominantPeriodWTSignificant, period] = significantCWT(x, dt, numRepeats);

periodPlot = hours(period);
figure(1013)
clf
subplot(1, 7, [1 2 3])
surf(t_WT/60, periodPlot, WT, 'EdgeColor', 'none')
title('Wavelet transform')
xlim([t_WT(1)/60 t_WT(end)/60])
ylim([periodPlot(1) periodPlot(end)])
load cmap
colormap(cmap)
colorbar
xlabel('Time (h)')
ylabel('Period (h)')
view(0,90)
set(gca,'yscale','log')

subplot(1, 7, [4 5 6])
surf(t_WT/60, periodPlot, significantWT, 'EdgeColor', 'none')
title('Wavelet transform with bootstrapping')
xlim([t_WT(1)/60 t_WT(end)/60])
ylim([periodPlot(1) periodPlot(end)])
colorbar
xlabel('Time (h)')
ylabel('Period (h)')

view(0,90)
set(gca,'yscale','log')

subplot(1, 7, 7)
semilogy(nanmean(significantWT, 2), periodPlot, 'k', 'linewidth', 1); hold on
ylim([periodPlot(1) periodPlot(end)])
drawnow

numRepeats = 1; %Number of times to test a random WT for significance
parfor cell = 1:cells
    x = P(cell, 0.2*Nt:end);
    [~, dominantPeriodWT(cell), ~, dominantPeriodWTSignificant(cell), period] = significantCWT(x, dt, numRepeats);
end

figure(1015)
clf
violinplot([hours(dominantPeriodWT)', hours(dominantPeriodWTSignificant)'], {'No significance test', 'Significance test'});

figure(16)
clf
meanExpression = mean(P(:, 0.2*Nt:end), 2);
scatter(dominantPeriodWTSignificant, meanExpression);
title('Wavelet period versus mean expression')
xlabel('Period (h)')
ylabel('Mean abundance')

end




