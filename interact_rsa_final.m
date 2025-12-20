function N = interact_rsa_final(x,p0,tmax)
%INTERACT_RSA  Interactive (polarization-dependent) 1D RSA simulation.
%
%   N = INTERACT_RSA(x,p0,tmax) runs a random sequential adsorption (RSA)
%   process on a 1D interval of initial length x. At each step, a unit-length
%   dipolar segment is attempted and deposited with orientation "+-" with
%   probability p0, or "-+" with probability (1-p0). The deposition rule
%   depends on the polarities at the boundaries of the interval where the
%   attempt occurs, which determines how the remaining free space is split.
%
%   Inputs
%   ------
%   x     : Size (length) of the initial interval (scalar).
%   p0    : Probability of adsorbing a dipole with orientation "+-" (0<=p0<=1).
%   tmax  : Maximum number of Monte Carlo steps (integer).
%
%   Output
%   ------
%   N     : Estimated total number of deposited dipoles at the end of the run.
%
%   Notes on representation
%   -----------------------
%   The system is stored in the cell array V as a concatenation of "interval
%   blocks" of length 3:
%       { left_polarity , free_length , right_polarity }
%   where left_polarity and right_polarity are '+' or '-', and free_length is
%   the available empty length inside that interval (a nonnegative scalar).
%
%   Each deposition attempt picks (implicitly) each current interval block in
%   sequence (ii = 1:ninterval) and, if there is enough free length (>=1), it
%   deposits one unit segment and replaces that interval by two new interval
%   blocks separated by the newly deposited dipole. The internal polarities
%   of the deposited dipole are inserted as either '+','-' or '-','+'.
%
%   Termination
%   -----------
%   The loop stops early if an iteration produces no net change in V
%   (i.e., length(V) stays constant), meaning no interval had free length >= 1,
%   so further adsorption is impossible.
%
%   Important
%   ---------
%   This function contains example parameter assignments (x, p0, tmax) inside
%   the function body. If you call the function with arguments, those example
%   assignments will override the input values. For public release, you may
%   want to remove or comment those example lines, but per your request the
%   structure is kept unchanged here.

% x: Size of the initial interval.
% p0: Probability of adsorv a dipole whit orientation +-.
% tmax: Max time (steps) of simultion

% Polarity of the initial interval is  +___-
% (i.e., left boundary '+', right boundary '-').

% Examples of values:
% NOTE: These lines override the function inputs if left active.
% x = 8;
% p0 = 1/2;
% tmax = 10000;

% V stores the current list of interval blocks. Each block has 3 entries:
%   {left_polarity, free_length, right_polarity}
% The initial condition is one interval: '+' --(length x)-- '-'.
V = {'+', x, '-'};

% ninterval counts how many interval blocks currently exist.
ninterval = 1;

% lengthV is used to detect when V stops changing (jamming reached).
lengthV = length(V);

for t = 1:tmax
    
    % 'plus' holds the starting indices (0-based offset) of each 3-entry block
    % within V. For interval ii, the corresponding entries are V{(1:3)+plus(ii)}.
    plus = (0:ninterval-1)*3;
    
    % newintervals will count how many interval blocks will exist after updating.
    newintervals = 0;
    
    for ii = 1:ninterval
         
        % Extract the ii-th interval block from the flattened cell array V.
        INTERVAL = {V{(1:3)+plus(ii)}};

        % If the free length is less than 1, no unit segment fits: keep interval.
        if INTERVAL{2}<1
            NEWINTERVAL{ii,:} = INTERVAL;
            newintervals = newintervals + 1;
            continue
        end

        % Decide the orientation of the next dipole:
        %   "+-" with probability p0, otherwise "-+".
        if rand<p0

            % Case: deposit "+-" (internal polarities '+','-').
            % The way the remaining free space is split depends on the boundary
            % polarities of the current interval (INTERVAL{1} and INTERVAL{3}).
            % INTERVAL{2} is the available free length in the interval.
            if strcmp(INTERVAL{1},'+') & strcmp(INTERVAL{3},'+') & INTERVAL{2}>=1
                % Same-sign boundaries (++): bias the remaining space to the left side.
                NEWINTERVAL{ii,:} = {INTERVAL{1}, INTERVAL{2}-1, '+', '-', 0 ,INTERVAL{3}};

            elseif strcmp(INTERVAL{1},'+') & strcmp(INTERVAL{3},'-') & INTERVAL{2}>=1
                % Opposite-sign boundaries (+-): split remaining space symmetrically.
                NEWINTERVAL{ii,:} = {INTERVAL{1}, (INTERVAL{2}-1)/2, '+', '-', (INTERVAL{2}-1)/2 ,INTERVAL{3}};

            elseif strcmp(INTERVAL{1},'-') & strcmp(INTERVAL{3},'+') & INTERVAL{2}>=1
                % Opposite-sign boundaries (-+): split remaining space symmetrically.
                NEWINTERVAL{ii,:} = {INTERVAL{1}, (INTERVAL{2}-1)/2, '+', '-', (INTERVAL{2}-1)/2 ,INTERVAL{3}};

            elseif strcmp(INTERVAL{1},'-') & strcmp(INTERVAL{3},'-') & INTERVAL{2}>=1
                % Same-sign boundaries (--): bias the remaining space to the right side.
                NEWINTERVAL{ii,:} = {INTERVAL{1}, 0, '+', '-', INTERVAL{2}-1, INTERVAL{3}};

            end

        else

            % Case: deposit "-+" (internal polarities '-','+').
            % Same logic as above, but with opposite dipole orientation.
            if strcmp(INTERVAL{1},'+') & strcmp(INTERVAL{3},'+') & INTERVAL{2}>=1
                NEWINTERVAL{ii,:} = {INTERVAL{1}, 0, '-','+', INTERVAL{2}-1 ,INTERVAL{3}};

            elseif strcmp(INTERVAL{1},'+') & strcmp(INTERVAL{3},'-') & INTERVAL{2}>=1
                NEWINTERVAL{ii,:} = {INTERVAL{1}, (INTERVAL{2}-1)/2, '-', '+', (INTERVAL{2}-1)/2 ,INTERVAL{3}};

            elseif strcmp(INTERVAL{1},'-') & strcmp(INTERVAL{3},'+') & INTERVAL{2}>=1
                NEWINTERVAL{ii,:} = {INTERVAL{1}, (INTERVAL{2}-1)/2, '-', '+', (INTERVAL{2}-1)/2 ,INTERVAL{3}};

            elseif strcmp(INTERVAL{1},'-') & strcmp(INTERVAL{3},'-') & INTERVAL{2}>=1
                NEWINTERVAL{ii,:} = {INTERVAL{1},  INTERVAL{2}-1, '-', '+', 0 , INTERVAL{3}};

            end

        end
        
        % If a deposition occurs, one interval becomes two intervals, so we add 2.
        newintervals = newintervals + 2;

    end

    % Rebuild V by concatenating all updated interval blocks stored in NEWINTERVAL.
    % This keeps the same flattened "triplet blocks" format.
    V = [];
    for ii = 1:ninterval
        V = [V NEWINTERVAL{ii,:}];
    end

    % If V did not change in total length, no adsorption took place anywhere:
    % the system is jammed, so stop.
    if lengthV == length(V)
        break
    end
    
    % Update counters for next iteration:
    ninterval = newintervals;     
    lengthV = length(V);

end

% N is returned as (newintervals - 1). In this encoding, newintervals is the
% number of interval blocks after the last processed step; subtracting 1 gives
% the number of deposited dipoles in this representation.
N = newintervals-1;
