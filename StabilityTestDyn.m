function [LMtest,cval,tmonth] = StabilityTestDyn(y,x,w0)

%% StabilityTestDyn.m
%-----------------------------------------------------------------------------------------------------------------------
%   This function performs the sup LM test and spits out the maximum value and the time at which the maximum value 
%   occurs. The function works for inclusion of up to 10 parameters in total and for the values of the trimming 
%   parameter given in Estrella (2003). 
%
%   Function inputs:
%   --------------------------------------------------------------------------------------------------------------------
%
%       y           = A T x 1 vector containing the binary dependent variable 
%       x           = A T x K matrix of explanatory variables. 
%       w0          = A scaler indicating how much to trim of the sample at each endpoint
%
%   Function outputs:
%   --------------------------------------------------------------------------------------------------------------------
%
%       LMtest      = Value of the test statistic
%       cval        = Critical value for the test statistic
%       tmonth      = Month in which the first break is observed
%
%   --------------------------------
%   Last modified: September 3, 2015
%   --------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------

% Checking for errors
if (nargin < 3)
    error('StabilityTestDyn.m: Not enough input parameters');
end

if (nargin > 3)
    error('StabilityTestDyn.m: Too many input parameters');
end

%-----------------------------------------------------------------------------------------------------------------------
%% GETTING FULL SAMPLE ESTIMATES
%-----------------------------------------------------------------------------------------------------------------------

% Get full sample parameter estimates
res     = dynProbit(y,x);
parm    = res.beta;

%-----------------------------------------------------------------------------------------------------------------------
%% CONSTRUCTING AUTOREGRESSIVE INDEX
%-----------------------------------------------------------------------------------------------------------------------

% Construct autoregressive index
x       = [ones(size(y,1),1) x];
t       = size(x,1);
b       = parm(1:end-1);
a       = parm(end);
pind    = nan(t,1);
cdf     = nan(t,1);
pdf     = nan(t,1);
pi0     = mean(x)*b / (1 - a);
pind(1) = pi0;
cdf(1)  = norm_cdf(pi0);
pdf(1)  = norm_pdf(pi0);
for iStep = 2:t
    pind(iStep) = a*pind(iStep-1) + x(iStep,:)*b;
    cdf(iStep)  = norm_cdf(pind(iStep));
    pdf(iStep)  = norm_pdf(pind(iStep));
end
h   = 1e-6; 
cdf = min(max(cdf,h),1-h);
x   = [x(2:end,:) pind(1:end-1)];

%-----------------------------------------------------------------------------------------------------------------------
%% COMPUTING ROBUST SCORE COVARIANCE MATRIX
%-----------------------------------------------------------------------------------------------------------------------

% Compute score function (gradiant)
yg      = y(2:end,1);
cdfg    = cdf(2:end);
pdfg    = pdf(2:end);
grad    = repmat((yg-cdfg)./(cdfg.*(1-cdfg)).*pdfg,1,length(parm)).*x;

% Compute robust score covariance matrix
[t, k]  = size(x);
nlag    = floor(4*(t/100)^(2/9));
d       = grad';
G       = zeros(k,k);
w       = zeros(2*nlag+1,1);
v       = 0;
while v ~= nlag+1

    ga  = zeros(k,k);
    z   = v / nlag;

    if z <= 0.5
     w(nlag+1+v,1) = 1 - 6*z^2 + 6*z^3;
    elseif z <= 1
         w(nlag+1+v,1) = 2*(1 - z)^3;
    else
         w(nlag+1+v,1) = 0;
    end

    za  = d(:,(v+1):t)*d(:,1:t-v)';

    if v == 0
      ga = ga + za;
    else
      ga = ga + za + za';
    end
    G   = G + w(nlag+1+v,1)*ga;
    v   = v + 1;

end

%-----------------------------------------------------------------------------------------------------------------------
%% COMPUTING THE LM TEST FOR BREAK POINTS
%-----------------------------------------------------------------------------------------------------------------------

% Compute LM over interior portion of sample
jj      = 0;
LMtest  = zeros(floor(t*(1-w0)) - floor(t*w0),1);
for iStep=floor(t*w0):1:floor(t*(1-w0))

    jj  = jj + 1;
    w1 = iStep/t;
    w2 = 1 - w1;
    LMtest(jj,1) = (1/(w1*w2)) * (sum(grad(1:iStep,:))/G)*sum(grad(1:iStep,:))';

end

% Compute test statistic and get critical value
[LMtest, tmonth]    = max(LMtest);
tmonth              = tmonth + floor(t*w0);
cval                = critval(w0,size(x,2));

end

%-----------------------------------------------------------------------------------------------------------------------
%% CRITICAL VALUES FROM ESTRELLA (2003) USING A 5% SIGNIFICANCE LEVEL
%-----------------------------------------------------------------------------------------------------------------------

function cval = critval(w0,p)

mcrit = [
    0.5 0.49 0.48 0.47 0.45 0.40 0.35 0.30 0.25 0.20 0.15 0.10 0.05
    3.84 4.73 5.09 5.37 5.80 6.57 7.14 7.61 8.04 8.45 8.86 9.31 9.90
    5.99 7.09 7.53 7.86 8.37 9.27 9.93 10.47 10.96 11.41 11.87 12.37 13.01
    7.81 9.06 9.55 9.92 10.49 11.49 12.21 12.80 13.33 13.82 14.31 14.85 15.53
    9.49 10.85 11.39 11.79 12.41 13.48 14.26 14.89 15.45 15.97 16.50 17.06 17.78
    11.07 12.54 13.12 13.55 14.21 15.35 16.17 16.84 17.42 17.98 18.52 19.12 19.87
    12.59 14.15 14.77 15.22 15.92 17.12 17.98 18.68 19.30 19.87 20.44 21.06 21.84
    14.07 15.71 16.36 16.84 17.57 18.83 19.72 20.45 21.09 21.69 22.28 22.92 23.73
    15.51 17.23 17.91 18.41 19.17 20.48 21.41 22.17 22.83 23.45 24.06 24.72 25.55
    16.92 18.72 19.42 19.94 20.74 22.09 23.06 23.84 24.52 25.16 25.79 26.47 27.31
    18.31 20.18 20.91 21.45 22.27 23.67 24.67 25.47 26.18 26.83 27.48 28.18 29.05
];
           
cval = mcrit(p+1,mcrit(1,:) == w0);
           
end

%-----------------------------------------------------------------------------------------------------------------------
% END OF FUNCTION
%-----------------------------------------------------------------------------------------------------------------------