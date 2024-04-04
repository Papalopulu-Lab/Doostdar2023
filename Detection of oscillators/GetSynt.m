function data=GetSynt(cov1,Noise,repeats,x)
CovMatrix1 = zeros(length(x),length(x));
for j = 1:length(x)
    CovMatrix1(:,j) = circshift(cov1,[j-1,0]);
end
CovMatrix1 = CovMatrix1';
CVM1 = triu(CovMatrix1,0)+(triu(CovMatrix1,1))';
MU = zeros(repeats,length(x));
Meas = diag((Noise^2).*ones(1,length(x)));
CVM1 = CVM1 + Meas;
SIGMA = CVM1; % change this to switch non-osc and osc
data= mvnrnd(MU,SIGMA);
%figure,plot(x,data(1,:))
end