function [S] = MC_HM(alpha, beta, D, periodFactRet)

% %%%%%%%Comment out test code once finished
% clc
% clear all
% format long
% data = load('matlab.mat');
% factorRet = getfield(data,'factorRet');
% periodFactRet = table2array(factorRet);
% D = csvread('D.csv');
% alpha = csvread('alpha.csv');
% beta = csvread('beta.csv');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
 %Calcualte correlation matrix
 corrMat = corrcov(Q_FF);
 
 %Correlation between factor column 1 and factor column 2
 u = copularnd('Gaussian',corrMat,5000);
 [s1, i1] =sort(u(:,1));
 [s2, i2] = sort(u(:,2));
 [s3, i3] = sort(u(:,3));
 
 x1 = zeros(size(s1));
 x2 = zeros(size(s2));
 x3 = zeros(size(s3));
 
 x1(i1) = sort(R(:,1));
 x2(i2) = sort(R(:,2));
 x3(i3) = sort(R(:,3));
 x = [x1 x2 x3];
S = alpha' + beta*x' + mvnrnd(zeros(20, 1), D, 5000)';
end