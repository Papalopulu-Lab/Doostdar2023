clear all, close all, clc
addpath(genpath('GPML'))
startup
addpath(genpath('raw data'))
%% Step 1. Choose lengthscale to detrend (in hours) and other parameters
Lengthscale =4.5; 
DetrendParam = log(1./(2*Lengthscale.^2));
%% HV Venus norm
fnames='190530_HV_VENUS_H2B_PBG.xls';
exptnames='190530_HV_VENUS_H2B_PBG';
coldata=[2:77]; 
bkgrdata=[78:80];
strow=4;  % start row is 4
%% HVP Venus norm
% fnames='190530_HVP_VENUS_H2B_PBG.xls';
% exptnames='190530_HVP_VENUS_H2B_PBG';
% coldata=[2:44]; 
% bkgrdata=[45:47];
% strow=4;  % start row is 4
%% HV keima nuclear marker
% fnames='190530_HV_keima_PBG.csv';
% exptnames='190530_HV_keima_PBG';
% coldata=[2:77]; 
% bkgrdata=[78:80];
% strow=4;  %start row is 4
%% HVP keima nuclear marker
% fnames='190530_HVP_keima_PBG.csv';
% exptnames='190530_HVP_keima_PBG';
% coldata=[2:44]; 
% bkgrdata=[45:47];
% strow=4;  % start row is 4
%% Step 2 -setup variables- check time is converted to h correctly 
% make an new results folder
dirname=[exptnames];
mkdir(dirname)
% read the data
num = xlsread(fnames);
num(isnan(num)) = 0;
data = num(strow:end,coldata);
bckgd= num(strow:end,bkgrdata);
time = num(strow:end,1); % enter column of time always 1
%
M1=zeros(1,size(data,2));
M2=zeros(1,size(data,2));
S1=zeros(1,size(data,2));
S2=zeros(1,size(data,2));
detrendData=[];
%% background SD
[stdev] = bckgdstdev(bckgd,time,-7); % this is around 0.1
%% Step 3. Fit individual cells
parfor i = 1:size(data,2)
    i
    y1 =data(:,i);
    x = time;
    x(y1==0) = []; %deletes times from which no signal
    samp = length(x);
    y1(y1==0) = [];
    raw = y1;
    % store normalization constants
    M1(i)=mean(y1);
    S1(i)=std(y1);
    %% attention here
    sd=stdev/std(y1);
    %y1 = (y1 - mean(y1))/std(y1);
    y1 = (y1 )/std(y1); % this was used why
    raw = y1;
    % detrend data 
    [m,m1,par0] = detrenddataNEW(raw,x,DetrendParam);
    y1 = y1-m; %detrended y1
    snr=var(raw)/sd;
    % save new normalization constants
    M2(i)=mean(y1);
    S2(i)=std(y1);
    y1 = (y1 )/std(y1);
    y1raw=y1*S1(i)*S2(i);
    detrendData=[detrendData; bring_to_size(y1raw',[1,numel(time)],NaN)];
    % fit OU and OUoscillatory models
    [LLR,LL1,LL2,par1,par2,xd,md]  = fitDataLLR(x,y1,1/snr);
    par0TOT(i,:) = par0; % trend
    par1TOT(i,:) = par1; % OU
    par2TOT(i,:) = par2; % OUosc
    LLRTOT(i) = LLR;
    LL1TOT(i)=LL1;
    LL2TOT(i)=LL2;
    dfit(i).x=x;
    dfit(i).xd=xd;
    dfit(i).Raw=raw; 
    dfit(i).RawOsc=y1; 
    dfit(i).M1=M1(i);
    dfit(i).M2=M2(i);
    dfit(i).S1=S1(i);
    dfit(i).S2=S2(i);
    dfit(i).M=m1*S1(i);
    dfit(i).m=md*S2(i)*S1(i)+M2(i);
    showfigure1(x,m,raw,y1,xd,md,M1(i),S1(i),S2(i),LLR,par2,i);
    print(gcf(),[dirname,'/Cell',num2str(i)],'-dpng');
end
%% plots distribution of LLR scores
figure()
histogram(LLRTOT,15) % choose number of bins with second number
xlabel('LLR score')
ylabel('Frequency')
title('Distribution of LLR scores')
print(gcf(),[dirname,'/LLR Scores'],'-dpng');
%% Hilbert Analysis with Fold change detection
hbt=HilbertFoldChange(dfit,dirname);
%%
close all
fc=[hbt(:).FoldChangeFit];
figure,boxplot(fc),title('Overall fold changes');
print(gcf(),[dirname,'/Overall fold change'],'-dpng');
for i=1:numel(hbt)
    try
        fcmax(i)=max([hbt(i).FoldChangeFit]);
    catch
        fcmax(i)=NaN;
    end
end
figure,boxplot(fcmax),title('Maximum fold change')
print(gcf(),[dirname,'/Maximum fold change'],'-dpng');
%% Step 4. Summary data plots 
% plot periods
periods = 2*pi()./par2TOT(:,2);
per=periods;
per(per>10)=[];
figure,histogram(per,5)
set(gca,'XLim',[0,3]);
title(['Periods ',num2str(round(mean(per)*100)/100),'h'])
xlabel('Period (hours)')
ylabel('Number of cells')
print(gcf(),[dirname,'/PeriodPlot'],'-dpng');
%% plots distribution of quality parameter of all cells
Quality=par2TOT(:,2)./(par2TOT(:,1));
% do not plot high qual
Quality(Quality>20)=[];
figure,histogram(Quality,5)
title('Quality of oscillations')
ylabel('Number of cells')
print(gcf(),[dirname,'/QualityPlot'],'-dpng');
%% Step 5. Export data analysis and save simulation files 
s=size(data);
savedata=-ones(s(2),9);
% read in the cell id
% Column 1 track id from file
savedata(:,1)=num(strow-1,coldata);
header{1}='TrackId';
% Column 2 processed cell id
savedata(:,2)=1:s(2);
header{2}='CellId';
% Column 3 lengthscale OU
savedata(:,3)=par1TOT(:,1);
header{3}='LenOU';
% Column 4 lengthscale OU
savedata(:,4)=par2TOT(:,1);
header{4}='LenOUosc';
% Column 5 period
savedata(:,5)=2*pi./par2TOT(:,2);
header{5}='Period(h)';
% Column 6 LLR
savedata(:,6)=LLRTOT;
header{6}='LLR';
% Column 7 quality
savedata(:,7)=(par2TOT(:,2)./(par2TOT(:,1)));
header{7}='Quality';
header{8}='FOD'; % !!saving done after globalFDR
% Column 8 fod label 1=osc 0=nonosc; initially set to -1
% savedata(:,8)=PassList;
header{9}='MaxFoldChange';
savedata(:,9)=fcmax;
if ismac==0
    xlswrite([dirname,'/SummaryFile'],header,'Sheet1','A1:I1');
    xlswrite([dirname,'/SummaryFile'],savedata,'Sheet2','A2');
else
    csvwrite([dirname,'/SummaryFile.csv'],savedata);
    writecell(header,[dirname, '/SummaryFileHeader.csv']);
    writecell(str(:),[dirname '/SummaryFileCellID.csv']);
end 
%% save simulation files here
save([dirname,'/',strtok(fnames,'.') ,'nofdr.mat'])
%% export complete Hilbert stats over time
FoldChangeMat=[];
PeakMat=[];
TroughMat=[];
PhaseMap=[];
for i=1:numel(hbt)
    vect=hbt(i).FoldChangeFit;
    FoldChangeMat=[FoldChangeMat,bring_to_size(vect',[100,1],NaN)];
    p=hbt(i).PairedIdx;
    PeakMat=[PeakMat, bring_to_size(p(:,1),[100,1],NaN)];
    TroughMat=[TroughMat,bring_to_size(p(:,2),[100,1],NaN)];
    ph=hbt(i).HilbertPhase;
    PhaseMap=[PhaseMap, bring_to_size(ph,[numel(time),1],NaN)];
end
if ismac==0
    xlswrite([dirname,'/Hilbert'],FoldChangeMat,'FoldChanges');
    xlswrite([dirname,'/Hilbert'],PeakMat,'HilbertPeaks');
    xlswrite([dirname,'/Hilbert'],TroughMat,'HilbertTroughs');
    xlswrite([dirname,'/Hilbert'],PhaseMap,'PhaseDiagram');
else
    csvwrite([dirname,'/Hilbert/FoldChanges.csv'],FoldChangeMat);
    csvwrite([dirname,'/Hilbert/HilbertPeaks.csv'],PeakMat);
    csvwrite([dirname,'/Hilbert/HilbertTroughs.csv'],TroughMat);
    csvwrite([dirname,'/Hilbert/PhaseDiagram.csv'],PhaseMap);
end



