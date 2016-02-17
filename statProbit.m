function result = statProbit(y,x)

%% statProbit.m
%-----------------------------------------------------------------------------------------------------------------------
%
%   This function estimates the parameters of a static probit model using maximum likelihood estimation. 
%   The function makes use of probitlikstat.m in the numerical optimization procedure. Standard errors are computed as 
%   suggested in Kauppi & Saikkonen (2008). An intercept is automatically included in the design matrix.
%
%   Function inputs:
%   --------------------------------------------------------------------------------------------------------------------
%
%       y           = A T x 1 vector containing the binary dependent variable 
%       x           = A T x K matrix of explanatory variables. 
%
%   Function outputs:
%   --------------------------------------------------------------------------------------------------------------------
%
%       result      = A structure of estimates and model diagnostics
%
%   --------------------------------
%   Last modified: September 3, 2015
%   --------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------

% Checking for correct number of arguments
if (nargin < 1)
    error('statProbit.m: Not enough input parameters'); 
end

if (nargin > 2) 
    error('statProbit.m: Too many input parameters'); 
end

% Checking for variation in the binary dependent variable
tmp     = find(y==1); 
chk     = length(tmp);
nobs    = size(y,1);
if (chk == nobs)
    error('statProbit.m: y-vector contains all ones');
elseif (chk == 0)
    error('statProbit.m: y-vector contains all zeros');
end

% Checking that y and x vector have same # of observations
nobss   = size(x,1);
if (nobs ~= nobss)
   error('statProbit.m: y- and x-vector are not equal length');
end

%-----------------------------------------------------------------------------------------------------------------------
%% MAXIMUM LIKELIHOOD ESTIMATION OF THE DYNAMIC PROBIT MODEL
%-----------------------------------------------------------------------------------------------------------------------

% Setting starting values for numerical optimization
iota    = ones(nobs,1);
x       = [iota x];
b       = (x'*x)\(x'*y);
[t, k]  = size(x);
parm    = b;

% Numerical optimization
opts = optimset('Display','off','MaxIter',3000,'MaxFunEvals',5000,'TolX',1e-6,'TolFun',1e-6,'LargeScale','off');
[parm,fval,eflag] = fminunc('probitlikstat',parm,opts,y,x);

%-----------------------------------------------------------------------------------------------------------------------
%% COMPUTING GRADIANT AND HESSIAN
%-----------------------------------------------------------------------------------------------------------------------

% Gradiant
cdf     = norm_cdf(x*parm);
h       = 1e-6; 
cdf     = min(max(cdf,h),1-h);
pdf     = norm_pdf(x*parm);
grad    = repmat((y-cdf)./(cdf.*(1-cdf)).*pdf, 1,length(parm)).*x;

% Hessian
q       = 2*y - iota;  
xparm   = x*parm;
pdfh    = norm_pdf(q.*xparm);
cdfh    = norm_cdf(q.*xparm);
lambda  = (q.*pdfh)./cdfh; 
H       = zeros(k,k);
for iStep = 1:t

    xb    = x(iStep,:)*parm;
    xp    = x(iStep,:)'; 
    H     = H - lambda(iStep,1)*(lambda(iStep,1) +xb)*(xp*x(iStep,:));

end

%-----------------------------------------------------------------------------------------------------------------------
%% PERFORMING PARZEN KERNEL CORRECTION (KAUPPI & SAIKKONEN, 2008)
%-----------------------------------------------------------------------------------------------------------------------

% Setting up preliminaries
nlag    = floor(4*(t/100)^(2/9));
d       = grad';
G       = zeros(k,k);
w       = zeros(2*nlag+1,1);
v       = 0;

% Computing covariance matrix
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

% Covariance estimators
%covparm = inv(grad'*grad); % BHHH Estimator
%covparm = inv(-H);         % The Hessian Estimator
covparm = H\G/H;            % The Kauppi & Saikkonen Estimator
stdparm = sqrt(diag(covparm));
tstat = parm ./ stdparm;

%-----------------------------------------------------------------------------------------------------------------------
%% FITTED RECESSION PROBABILITIES AND OTHER USEFULL STATISTICS
%-----------------------------------------------------------------------------------------------------------------------

% Fitted recession probabilities
yhat    = cdf;

% Residuals
e       = y - yhat;
 
% Computing sample recessions proportions
temp    = find(y==1);
P       = length(temp);
P       = P/nobs;

% Computing the log likelihood ratio
like0   = nobs*(P*log(P) + (1-P)*log(1-P)); % Restricted likelihood
like1   = -fval; % Unrestricted likelihood
llratio = 2*(like1-like0);

% Computing information criteria
AIC     = -(like1 - k);
BIC     = -(like1 - 0.5*k*log(nobs));

%Computing pseudo-R2 measures
r2mcf   =  1 - (abs(like1)/abs(like0));         %McFadden
r2est   = 1-(like1/like0)^(-(2/nobs)*like0);    %Estrella

% Computing forecast evaluation measures
QPS     = (2/nobs) * sum((yhat-y).^2);
LPS     = (-1/nobs) * sum( y .* log(yhat) + (1-y) .* log(1-yhat));
KS      = (sum(y.*(yhat>0.5)) / (sum(y))) - (sum((1-y).*(yhat>0.5)) / sum(1-y));

% Computing proportion of correct classifications
indx    = find(y == 1);
CR50    = sum(((yhat(indx,1) > 0.5) == y(indx,1)))/size(y(indx,1),1);
CR25    = sum(((yhat(indx,1) > 0.25) == y(indx,1)))/size(y(indx,1),1);
indx    = find(y == 0);
CE50    = sum(((yhat(indx,1) > 0.5) == y(indx,1)))/size(y(indx,1),1);
CE25    = sum(((yhat(indx,1) > 0.25) == y(indx,1)))/size(y(indx,1),1);

% Computing AUROC
[Xa,Ya,Ta,AUC]  = perfcurve(y,yhat,1);
Q1              = AUC/(2-AUC);
Q2              = (2*AUC^2) / ( 1 + AUC);
n0              = size(y(indx,1),1);
n1              = nobs-n0;
AUCsigma2       = (1/(n0*n1))*(AUC*(1-AUC)+(n1-1)*(Q1-AUC^2)+(n0-1)*(Q2-AUC^2));
AUCsigma        = sqrt(AUCsigma2);
AUCtest         = (AUC-0.5)/AUCsigma;
AUCpval         = (1-normcdf(AUCtest,0,1));

% Computing optimal cut-off point
[~,maxbj2]  = max((2.*mean(y).*Ya-mean(y)) - (2.*(1-mean(y)).*Xa - (1-mean(y))));
optroc      = Ta(maxbj2);

% Saving the important results
result.beta     = parm;
result.resid    = e;
result.yhat     = yhat;
result.stdb     = stdparm;
result.tstat    = tstat;
result.bic      = BIC;
result.aic      = AIC;
result.likhood  = like1;
result.llratio  = llratio;
result.r2mcf    = r2mcf;
result.r2est    = r2est;
result.nber     = y;
result.qps      = QPS;
result.lps      = LPS;
result.ks       = KS;
result.auc      = AUC;
result.aucse    = AUCsigma;
result.auctest  = AUCtest;
result.aucp     = AUCpval;
result.optroc   = optroc;
result.cr50     = CR50;
result.cr25     = CR25;
result.ce50     = CE50;
result.ce25     = CE25;
result.flag     = eflag;

%-----------------------------------------------------------------------------------------------------------------------
%% CALCULATING THE NORMAL PROBABILITY DENSITY FUNCTION
%-----------------------------------------------------------------------------------------------------------------------

function pdf = norm_pdf(x,m,v)

% Error checking
if ~(nargin == 1) || (nargin == 3)
    error('Wrong # of arguments to norm_pdf');
elseif (nargin == 1) % standard normal with mean 0 and variance 1
    pdf         = stdn_pdf(x);
else
    r           = size(x,1);
    pdf         = zeros(r,1);
    pdf(1:r,1)  = stdn_pdf((x(1:r,1) - m(1:r,1)) ./ sqrt(v(1:r,1))) ./ sqrt(v(1:r,1));
end

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
%% CALCULATING THE STANDARD NORMAL PROBABILITY DENSITY FUNCTION
%-----------------------------------------------------------------------------------------------------------------------

function pdf = stdn_pdf(x)

% Error checking
if (nargin ~= 1)
    error('Wrong # of arguments to stdn_pdf');
end

[r, c]  = size(x);
s       = r * c;
x       = reshape(x,1,s);
pdf     = zeros(1,s);

k = find(isnan(x));
if any(k)
    pdf(k) = NaN * ones(1,length(k));
end

k = find(~isinf(x));
if any(k)
    pdf(k) = (2 * pi)^(-1/2) * exp(-x(k).^2 / 2);
end

pdf = reshape(pdf,r,c);

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