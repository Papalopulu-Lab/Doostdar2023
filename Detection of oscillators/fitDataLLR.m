function [LLR,LL1,LL2,par1,par2,xt,m] = fitDataLLR(x,y1,a0)
N=10; % number of runs
repeats=N;
LLR = zeros(N,1);
par1vec = zeros(N,3);
par2vec = zeros(N,4);
% init lengthscale keep only positive values
avec=0.5+0.05*randn(100,1); avec(avec<0)=[]; avec=avec(1:N); % 0.5 works ok
%% optimise OU model Kou=sigamsig*exp(-avec*t)+sigmanoise
for i = 1:N;
    samp = length(x);
    likfunc = @likGauss;
    par1mean = [log(avec(i)),log(var(y1))]; % was 0.5
    % settings
    covfunc = @covOUa;
    hyp1.lik =log(a0);
    hyp1.cov = par1mean;
    %prior on noise estimated to be 10-15% of signal mean
    %prior.cov={{@priorSmoothBox1,log(0.5),log(2),10};[];}; % prior for lengthscale
    prior.lik ={{@priorDelta}};
    inf = {@infPrior,@infExact,[]};
    hyp1 = minimize(hyp1, @gp, -1000, inf, [], covfunc, likfunc, x, y1);
    % negative log likelihood- invert sign for log likelihood
    nlmlOU = gp(hyp1, @infExact, [], covfunc, likfunc, x, y1);
    % matches mean of large number of realisations-no crossval needed here
    % statistic used is 2log(ll)
    LLou(i)=-2*nlmlOU/numel(x);
    % order of estimated parameters: alpha, sigmasig sigmanoise
    hypOU(i)=hyp1;
    clear prior
end
% get best OU model
[LL1,idx1] =max(LLou);
par1 = [exp(hypOU(idx1).cov),exp(hypOU(idx1).lik)];
%% fit the OU osc model with fixed noise
% start solution around the same cond - use avec from above
%avec = 0.1 + (0.2)*randn(100,1); avec(avec<0)=[]; avec=avec(1:N);
bvec = 5 + (0.5)*randn(N,1); %2 + (4)*rand(N,1); % period
clear hyp2
for i=1:N
    % use estimates from this model to init the next
    covfunc = @covOUosca; 
    par2mean = [log(avec(i)),log(2*pi/bvec(i)),log(var(y1))];
    hyp2.cov = par2mean;
    prior.cov={{@priorSmoothBox1,log(1),log(2),5};[];[];}; % weak prior on lengthscale
    hyp2.lik=log(a0);
    prior.lik={{@priorDelta}};
    inf = {@infPrior,@infExact,[]};
    hyp2 = minimize(hyp2, @gp, -1000, inf, [], covfunc, likfunc, x, y1); 
    nlmlOSC = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1);
    LLosc(i)=-2*nlmlOSC/numel(x); % checked to see if they match cross
    hypOUosc(i)=hyp2;
    clear prior
end
% choose best OU osc
[LL2,idx2] =max(LLosc);
par2 = [exp(hypOUosc(idx2).cov),exp(hypOUosc(idx2).lik)];
% generate fitting here
xt=[x(1):(x(2)-x(1))/20:x(end)]';
[m,s2] = gp(hypOUosc(idx2),inf,[],@covOUosca,@likGauss,x,y1,xt);
%% Log-likelihood ratio 
if LL2>LL1
    LLR=(LL2-LL1)*100;
else
    LLR=0;
end
end
