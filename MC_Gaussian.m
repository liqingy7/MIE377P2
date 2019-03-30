function [S] = MC_Gaussian(alpha, beta, D, periodFactRet)

    mu_FF = geomean(periodFactRet+1)-1;
    Q_FF = cov(periodFactRet);
% %This function calculates stock price scenarios using the Gaussian Monte
% %Carlo simulations
    S =  alpha' + beta*mvnrnd(mu_FF, Q_FF, 5000)' + mvnrnd(zeros(20, 1), D, 5000)';

end