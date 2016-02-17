%----------------------------------------------------------------------------------------------
% CALCULATING THE NORMAL PROBABILITY DENSITY FUNCTION
%----------------------------------------------------------------------------------------------

function pdf = norm_pdf(x,m,v)

% Error checking
if ~(nargin == 1) || (nargin == 3)
    error('Wrong # of arguments to norm_pdf');
elseif (nargin == 1) % standard normal with mean 0 and variance 1
    pdf = stdn_pdf(x);
else
    r = size(x,1);
    pdf = zeros(r,1);
    pdf(1:r,1) = stdn_pdf((x(1:r,1) - m(1:r,1)) ./ sqrt(v(1:r,1))) ./ sqrt(v(1:r,1));
end