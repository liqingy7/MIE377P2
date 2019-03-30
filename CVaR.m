function [x] = CVaR(rets)
    % Defining parameters
    l = 0.1;
    alpha = 0.95;
    S = 5000; 
    n = 20; 
    
    mu = ( geomean(rets + 1) - 1 )';     
    
    % Defining bounds
    lb = [];
    ub = [];

    
    % Define the inequality constraint matrices A and b
    A = -[zeros(S, n) eye(S) zeros(S, 1);rets eye(S) ones(S,1)];
    b = [zeros(S, 1); zeros(S, 1)];

    % Define the equality constraint matrices A_eq and b_eq
    Aeq = [ones(1,n) zeros(1,S) 0];
    beq = 1;

    % Define our objective linear cost function c
    k = (1 / ((1 - alpha) * S));
    c =  [-l.*mu;k*ones(S,1);1];


    % Use 'linprog' to find the optimal portfolio
    y = linprog( c, A, b, Aeq, beq, lb, ub );

    % Retrieve the optimal portfolio weights
    x = y(1:n);

end

