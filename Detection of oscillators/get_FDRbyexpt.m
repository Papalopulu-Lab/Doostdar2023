clear all, close all, clc
% put the data together and generate synth data to get a cutoff
LLRmat=[];
par1mat=[];
%% 190228 HV 
load('190530_HV_VENUS_H2B_PBG_vbedited_4.5_20-Oct-2022\190530_HV_VENUS_H2B_PBG_vbeditednofdr.mat','par1TOT','LLRTOT','time','dirname');
par1mat=[par1mat; par1TOT];
LLRmat=[LLRmat;LLRTOT'];
n1=numel(LLRTOT);
dirname1=dirname;
%time=time(1:numel(time));
%%
load('190530_HV_keima_PBG_vbedit_4.5_20-Oct-2022\190530_HV_keima_PBG_vbeditnofdr.mat','par1TOT','LLRTOT','time','dirname');
par1mat=[par1mat; par1TOT];
LLRmat=[LLRmat; LLRTOT'];
dirname2=dirname;
%%
load('190530_HV_venus_PBG_vbedited_4.5_20-Oct-2022\190530_HV_venus_PBG_vbeditednofdr.mat','par1TOT','LLRTOT','time','dirname');
par1mat=[par1mat; par1TOT];
LLRmat=[LLRmat; LLRTOT'];
dirname2=dirname;
%% generate synthetic OU models
par1M=par1mat;
repeats=20;%20 repeats per cell
dataTOT=[];
par1TOT=[];
for i = 1:size(par1M,1) % for each cell
    cov1=par1M(i,2)*exp(-par1M(i,1)*time);
    Noise=par1M(i,3);
    data=GetSynt(cov1,Noise,repeats,time);
    dataTOT=[dataTOT;data];
    par1TOT=[par1TOT;repmat(par1M(i,:),repeats,1)];
end
%%
parfor i=1:size(dataTOT,1)
    y1=dataTOT(i,:)';
    x=time;
    snr=var(y1)/0.1;
    a0=1/snr;
    [LLRs,LL1s,LL2s,par1s,par2s] = fitDataLLR(x,y1,a0);
    par1M(i,:) = par1s;
    par2M(i,:) = par2s;
    LLRS(i,:) = LLRs;
end
%% plot the distributions
LLRM=LLRmat;
grid=[0:2:max(LLRM)];
figure,subplot(2,1,1)
histogram(LLRM,grid), title('Her6 data')
subplot(2,1,2)
histogram(LLRS,grid),title('Synthetic data')
%% FDR part here
BICdiffTOT=LLRM;
BICdiffsynthTOT=LLRS;
upper = max([BICdiffTOT;BICdiffsynthTOT]);
lower1 = min([BICdiffTOT;BICdiffsynthTOT]);
lower = upper - 0.9*(upper-lower1); % 0.75 before usually 0.9
range = linspace(lower,upper,20);
for i= 1:length(range)
%     i
    cutoff = range(i);
    num = sum(BICdiffTOT<cutoff)/length(BICdiffTOT); % number of nonosc cells
    denom = sum(BICdiffsynthTOT<cutoff)/length(BICdiffsynthTOT); % number of nonosc cells in synthetic OU data
    piest(i) =  num/denom;
end
figure()
subplot(1,2,1)
plot(range,piest)
%ylim([0 1])
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
%
% idx=find(piest==0);
% piest(idx)=[];
% range(idx)=[];
%%cubic spline regression
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
subplot(1,2,2),plot(xx,yy)
%plot piGuess
hold on
plot(xx,yy,'color','r')
%ylim([0 1])
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
a = xlim();
b = ylim();
% text(a(1)-0.1*(a(2)-a(1)),b(2)+0.05*(b(2)-b(1)),{'a)'},...
%     'FontSize',9,'HorizontalAlignment','right', 'VerticalAlignment','bottom','color','k')
% do not accept a value of zero
piGUESS1= min(1,min(yy)); % was min(yy)
%% Go through BIClist calculating q values
binnings = 0:0.2:max([BICdiffTOT;BICdiffsynthTOT]);
q1 = zeros(length(binnings),1);
%find BIC list location above thresh...
for i = 1:length(binnings)
    Thresh = binnings(i);
%     (sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT));
%     (sum(BICdiffTOTs>Thresh)/length(BICdiffTOTs));
    q1(i) = piGUESS1*(sum(BICdiffsynthTOT>=Thresh)/length(BICdiffsynthTOT))/(sum(BICdiffTOT>=Thresh)/length(BICdiffTOT));
end
q = 0.03; % choose initial cutoff level
cutoff = find(q1<=q,1,'first')
Thresh=binnings(cutoff)
%
PassList=zeros(size(BICdiffTOT));
try
    PassList = BICdiffTOT>=Thresh;
end
%%
save globalFDR_190530HV_all.mat
%% save to summary files
% read the summary stats
savedata=xlsread([dirname1,'/SummaryFile.xls']);
savedata(:,8)=PassList(1:n1);
xlswrite([dirname1,'/SummaryFileFDR.xls'],savedata);
% 
savedata=csvread([dirname2,'/SummaryFile.xls']);
savedata(:,8)=PassList(n1+1:end);
xlswrite([dirname1,'/SummaryFileFDR.xls'],savedata);