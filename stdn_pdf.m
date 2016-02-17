%----------------------------------------------------------------------------------------------
% CALCULATING THE STANDARD NORMAL PROBABILITY DENSITY FUNCTION
%----------------------------------------------------------------------------------------------

function pdf = stdn_pdf(x)

% Error checking
if (nargin ~= 1)
    error('Wrong # of arguments to stdn_pdf');
end

[r, c] = size(x);
s = r * c;
x = reshape(x,1,s);
pdf = zeros(1,s);

k = find(isnan(x));
if any(k)
    pdf(k) = NaN * ones(1,length(k));
end

k = find(~isinf(x));
if any(k)
    pdf(k) = (2 * pi)^(-1/2) * exp(-x(k).^2 / 2);
end

pdf = reshape(pdf,r,c);