clear all, close all, clc
% put the HV and HVP data together and generate synth data to get a cutoff
LLRmat=[];
par1mat=[];
%% 190228 HV 
load('190530_HV_VENUS_H2B_PBG_27-Mar-2024\190530_HV_VENUS_H2B_PBGnofdr.mat','par1TOT','LLRTOT','time','dirname');
par1mat=[par1mat; par1TOT];
LLRmat=[LLRmat;LLRTOT'];
n1=numel(LLRTOT);
dirname1=dirname;
%%
load('190530_HVP_VENUS_H2B_PBG_27-Mar-2024\190530_HVP_VENUS_H2B_PBGnofdr.mat','par1TOT','LLRTOT','time','dirname');
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
upper = max([LLRM;LLRS]);
lower1 = min([LLRM;LLRS]);
lower = upper - 0.85*(upper-lower1); 
range = linspace(lower,upper,20);
for i= 1:length(range)
%     i
    cutoff = range(i);
    num = sum(LLRM<cutoff)/length(LLRM); % number of nonosc cells
    denom = sum(LLRS<cutoff)/length(LLRS); % number of nonosc cells in synthetic OU data
    piest(i) =  num/denom;
end
figure()
subplot(1,2,1)
plot(range,piest)
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
xx = linspace(lower,upper,100);
yy = spline(range,piest,xx);
subplot(1,2,2),plot(xx,yy)
hold on
plot(xx,yy,'color','r')
hold off
xlabel('\lambda')
ylabel('\pi_0(\lambda)')
xlim([lower1,upper])
a = xlim();
b = ylim();
% do not accept a value of zero
piGUESS1= min(1,min(yy)); 
%% Calculate q values
binnings = 0:0.2:max([LLRM;LLRS]);
q1 = zeros(length(binnings),1);
for i = 1:length(binnings)
    Thresh = binnings(i);
    q1(i) = piGUESS1*(sum(LLLRS>=Thresh)/length(LLRS))/(sum(LLRM>=Thresh)/length(LLRM));
end
q = 0.03; % choose initial cutoff level
cutoff = find(q1<=q,1,'first')
Thresh=binnings(cutoff)
%
PassList=zeros(size(LLRM));
try
    PassList = LLRM>=Thresh;
end
disp ('% Oscillators in dataset 1 -HV');
sum(PassList(1:n1))/numel(PassList(1:n1))
disp('% Oscillators in dataset 2- HVP')
sum(PassList(n1+1:end))/numel(PassList(n1+1:end))
%%
save globalFDR_190530HV_HVP_VenusH2B.mat
%% save to summary files
[tmp,header]=xlsread([dirname1,'/SummaryFile.xls'],'Sheet1');
xlswrite([dirname1,'/SummaryFileFDR.xls'],header,'Sheet1')
savedata=xlsread([dirname1,'/SummaryFile.xls'],'Sheet2');
savedata(:,8)=PassList(1:n1);
xlswrite([dirname1,'/SummaryFileFDR.xls'],savedata,'Sheet2','A2');
% 
[tmp,header]=xlsread([dirname2,'/SummaryFile.xls'],'Sheet1');
xlswrite([dirname2,'/SummaryFileFDR.xls'],header,'Sheet1');
savedata=xlsread([dirname2,'/SummaryFile.xls'],'Sheet2');
savedata(:,8)=PassList(n1+1:end);
xlswrite([dirname2,'/SummaryFileFDR.xls'],savedata,'Sheet2','A2');