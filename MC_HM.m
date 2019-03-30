%function [S] = MC_HM(alpha, beta, D, periodFactRet)

%%%%%%%Comment out test code once finished
clc
clear all
format long
data = load('matlab.mat');
factorRet = getfield(data,'factorRet');
periodFactRet = table2array(factorRet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 skew = skewness(periodFactRet);
 kurt = kurtosis(periodFactRet);
 mu_FF = geomean(periodFactRet+1)-1;
    Q_FF = cov(periodFactRet);
 %q = mvnrnd(mu_FF, Q_FF, 5000)'
 sigma = sqrt(diag(Q_FF));
R = zeros(5000,size(periodFactRet,2));
 for i = 1:size(periodFactRet,2)   
    R(:,i) = pearsrnd(mu_FF(i),sigma(i),skew(i),kurt(i),5000,1);
 end
%end