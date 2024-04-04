function [m,m1,par0] = detrenddataNEW(y1,x,lscale);
%DETRENDDATA Summary of this function goes here
%   Detailed explanation goes here
%%
likfunc = @likGauss; 
covfunc = @covSEa;
hyp2.cov = [lscale-1,log(var(y1))];
hyp2.lik = log(var(y1));

x1 = -10; x2 = lscale; eta = 16;
prior.cov = {{@priorSmoothBox2,x1,x2,eta}; {}};
inf = {@infPrior,@infExact,prior};
% inf = {@infExact,@infPrior,prior};

% hyp2 = minimize(hyp2, @gp, -10000, @infExact, [], covfunc, likfunc, x, y1);
hyp2 = minimize(hyp2, @gp, -10000, inf, [], covfunc, likfunc, x, y1);
par0 = [exp(hyp2.cov),exp(hyp2.lik)];

%%
z = x;%linspace(0, x(end), 1000)';
xt=[x(1):(x(2)-x(1))/20:x(end)]';
[m s2] = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1, z); % was z
[m1 s12] = gp(hyp2, @infExact, [], covfunc, likfunc, x, y1, xt); % was z
clear prior
%%
end

