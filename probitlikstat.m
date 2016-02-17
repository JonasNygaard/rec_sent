function like = probitlikstat(parm,y,x)

%% probitlikstat.m
%-----------------------------------------------------------------------------------------------------------------------
%   This function computes the value of the log likelihood function for a standard static probit model.
%
%   Function inputs:
%   --------------------------------------------------------------------------------------------------------------------
%
%       parm        = Vector of parameters of the model to be estimated
%       y           = A T x 1 vector containing the binary dependent variable 
%       x           = A T x K matrix of explanatory variables. 
%
%   Function outputs:
%   --------------------------------------------------------------------------------------------------------------------
%
%       like        = Value of the log likelihood function
%
%   --------------------------------
%   Last modified: September 3, 2015
%   --------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------

% Error checking
if (nargin < 3)
    error('probitlikstat.m: Not enough input parameters');
end

if (nargin > 3)
    error('probitlikstat.m: Too many input parameters');
end

smt = size(parm,2);
if smt ~= 1
    error('probitlikstat.m: Requires parm to be column vector');
end

%-----------------------------------------------------------------------------------------------------------------------
%% PRELIMINARY SETTINGS
%-----------------------------------------------------------------------------------------------------------------------

% Setting preliminaries
tol     = 1e-6;
iot     = ones(length(y),1);
cdf     = norm_cdf(x*parm);

%-----------------------------------------------------------------------------------------------------------------------
%% COMPUTING LOG LIKELIHOOD
%-----------------------------------------------------------------------------------------------------------------------

% Computing value of log likelihood
tmp     = find(cdf <= 0);
n1      = size(tmp,1);
if n1 ~= 0
    cdf(tmp,1) = tol*ones(length(tmp),1);
end

tmp     = find(cdf >= 1);
n1      = size(tmp,1);
if n1 ~= 0
    cdf(tmp,1) = (1-tol)*ones(length(tmp),1);
end

out     = y.*log(cdf)+(iot-y).*log(iot-cdf);
like    = -sum(out);

%-----------------------------------------------------------------------------------------------------------------------
%% CALCULATING THE NORMAL CUMULATIVE DISTRIBUTION FUNCTION
%-----------------------------------------------------------------------------------------------------------------------

function cdf = norm_cdf(x,m,v)

% Error checking
[r, c] = size(x);
if (r*c == 0)
    error('x must not be empty in norm_cdf');
end

if (nargin ==1)
    m = zeros(r,1);
    v = ones(r,1);
end

cdf = zeros(r,1);
cdf(1:r,1) = stdn_cdf((x(1:r,1) - m(1:r,1)) ./ sqrt(v(1:r,1)));

%-----------------------------------------------------------------------------------------------------------------------
%% CALCULATING THE STANDARD NORMAL CUMULATIVE DISTRIBUTION FUNCTION
%-----------------------------------------------------------------------------------------------------------------------

function cdf = stdn_cdf(x)

% Error checking
if (nargin ~= 1)
    error('Wrong # of arguments to stdn_cdf');
end

cdf = 0.5*(ones(size(x))+erf(x/sqrt(2)));

%-----------------------------------------------------------------------------------------------------------------------
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------------------------------