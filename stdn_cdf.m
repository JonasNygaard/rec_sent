%----------------------------------------------------------------------------------------------
% CALCULATING THE STANDARD NORMAL CUMULATIVE DISTRIBUTION FUNCTION
%----------------------------------------------------------------------------------------------

function cdf = stdn_cdf(x)

% Error checking
if (nargin ~= 1)
    error('Wrong # of arguments to stdn_cdf');
end

cdf = 0.5*(ones(size(x))+erf(x/sqrt(2)));