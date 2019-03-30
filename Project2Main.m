%% MIE377 (Winter 2019) - Project 1
% The purpose of this program is to implementtwo different stochastic 
% processes to simulate scenarios: 
% (i) a Gaussian process, and 
% (ii) a non-normal stochastic process with higher moments
% 
% and to use these simulations to test the out-of-sample performance of two
% portfolio optimization models: 
% 1. CVaR
% 2. Robust CVaR
%
% Student Name: Joshua Chang, Jinansh Shah 
% Student ID: 1003083147, 1003062614



clc
clear all
format short

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



adjClose = readtable('Project1_Data_adjClose.csv');
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
adjClose.Date = [];

% Load the factors weekly returns
factorRet = readtable('Project1_Data_FF_factors.csv');
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Date));
factorRet.Date = [];

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial budget to invest
initialVal = 100;

% Start of in-sample calibration period 
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Monte-Carlo simulations 
%MCList = {'MC_Gaussian' 'MC_HM'};
MCList = {'MC_Gaussian'};
MCList = cellfun(@str2func, MCList, 'UniformOutput', false);
NoSimulations = length(MCList);


% Investment strategies
invList = {'CVaR' 'CVaR_robust'};
invList = cellfun(@str2func, invList, 'UniformOutput', false);
NoStrats = length(invList);

n = 20; % number of assets
S = 5000; % number of scenarios


% Tags for factor models under different investment strategies
% tags = {'MVO (CAPM)' 'Card MVO (CAPM)' 'MVO (FF)' 'Card MVO (FF)' ...
%         'MVO (PCA)' 'Card MVO (PCA)'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period,
% transaction cost variable, and current value of the portfolio.
currentVal = zeros(NoPeriods, NoStrats * NoSimulations); 
toDay = 0;

for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        currentVal= ones(1, 6);
        currentVal(t,:) = initialVal;
        
        
    else
        for k = 1 : (NoStrats * NoSimulations)
            temp = transpose(NoShares{k})* currentPrices;
            currentVal(t,k) = temp;
            
        end
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    [alpha, beta, D, mu, Q] = FF(periodReturns, periodFactRet);
    
    for i = 1:NoSimulations
        rets{i} = MCList{i}(alpha, beta, D, periodFactRet); 
        for j = 1:NoStrats
            k = j + (i - 1) * 2;
            x{k}(:,t) = invList{j}(rets{i}');
        end
    end
    
    
    % Calculate the optimal number of shares of each stock you should hold
    for k = 1 : NoStrats * NoSimulations
        % Number of shares your portfolio holds per stock
        NoShares{k} = x{k}(:,t) .* currentVal(t,k) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,k) = periodPrices * NoShares{k};
        
        NoSharesOld{k} = NoShares{k};
        %------------------------------------------------------------------
        
    end

    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);

end

