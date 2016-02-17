%----------------------------------------------------------------------------------------------
% CALCULATING THE NORMAL CUMULATIVE DISTRIBUTION FUNCTION
%----------------------------------------------------------------------------------------------

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