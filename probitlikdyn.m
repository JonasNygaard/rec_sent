function like = probitlikdyn(parm,y,x)

%% probitlikdyn.m
%-----------------------------------------------------------------------------------------------------------------------
%   This function computes the value of the log likelihood function for a dynamic probit model that contains an 
%   autoregressive component.
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
    error('probitlikdyn.m: Not enough input parameters');
end

if (nargin > 3)
    error('probitlikdyn.m: Too many input parameters');
end

smt = size(parm,2);
if smt ~= 1
    error('probitlikdyn.m: Requires parm to be column vector');
end

%-----------------------------------------------------------------------------------------------------------------------
%% PRELIMINARY SETTINGS
%-----------------------------------------------------------------------------------------------------------------------

% Setting preliminaries
tol = 1e-6;
t   = size(y,1);
b   = parm(1:end-1);
a   = parm(end);

%-----------------------------------------------------------------------------------------------------------------------
%% CONSTRUCTING AUTOREGRESSIVE INDEX
%-----------------------------------------------------------------------------------------------------------------------

% Preallocations
pind    = nan(t,1);
cdf     = nan(t,1);

% Constructing the autoregressive index
pi0     = mean(x)*b / (1-a);
pind(1) = pi0;
cdf(1)  = norm_cdf(pi0);
for iStep = 2:t
    pind(iStep) = a*pind(iStep-1) + x(iStep,:)*b;
    cdf(iStep)  = norm_cdf(pind(iStep));
end
cdf     = min(max(cdf, tol),1-tol); 
out     = y.*log(cdf)+(1-y).*log(1-cdf);
out     = out(2:t);
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