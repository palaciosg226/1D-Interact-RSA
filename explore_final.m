clc
clear all

% -------------------------------------------------------------------------
% Batch driver script: run many independent RSA replicas for multiple x values
% -------------------------------------------------------------------------
% This script scans the initial interval length x over the vector X and, for
% each x, runs Nrep independent Monte Carlo replicas using interact_rsa().
% The outputs are aggregated into mean, variance, and a normalized standard
% deviation, then saved to a .mat file after each x (useful for long runs).
% -------------------------------------------------------------------------

% Probability of depositing a dipole with orientation "+-"
p0 = 1/2;

% Complementary probability for "-+" orientation (not used explicitly below,
% but kept here for clarity/documentation)
q0 = 1-p0;

% Range of initial interval lengths (scanned parameter)
% NOTE: This includes x = 0; depending on interact_rsa() implementation,
% x=0 may be trivial and/or can affect normalization (see STD_THETA below).
X = 0.02:0.02:63;

% Number of independent replicas (Monte Carlo runs) for each x
Nrep = 2e6;

% Maximum number of Monte Carlo steps per replica
tmax = 100000;

% Counter used only for progress reporting
c = 0;

for x = X
    c = c + 1;

    % Clear command window and show progress percentage
    clc
    disp([num2str(c/length(X)*100) ' % ...'])
    
    % Preallocation is recommended for performance, but structure is kept
    % unchanged as requested. The vector 'n' stores the output N from each replica.
    
    % Parallel loop over replicas:
    % Requires Parallel Computing Toolbox and an active parallel pool.
    % Each iteration calls interact_rsa(x,p0,tmax) and stores the result.
    parfor rep = 1:Nrep
       n(rep) = interact_rsa(x,p0,tmax);
    end

    % Mean number of deposited dipoles for this x
    N(c) = mean(n);

    % Standard deviation of the coverage fraction
    STD_THETA(c) = std(n/x);

    % Variance of the number of deposited dipoles for this x
    VAR_N(c) = var(n);

    % Save partial results after each x so you don't lose progress if the job stops.
    % This overwrites the same file each time; it will contain the latest arrays.
    save('results.mat')
end
