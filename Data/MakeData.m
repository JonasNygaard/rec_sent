%% MakeData.m
%-----------------------------------------------------------------------------------------------------------------------
%{
    This script loads and generates the data for the paper "Forecasting US Recessions: The Role of Sentiment". 

    --------------------------------
    Last modified: December, 2015
    --------------------------------
%}
%-----------------------------------------------------------------------------------------------------------------------

clear; clc; close all; tStart = tic; format shortg; c = clock; addpath('../');
disp('-----------------------------------------------------------------------------');
disp('Running the MakeData.m script.                                               ');
fprintf('Code initiated at %.0f:%.0f on %.0f / %0.f - %0.f \n',c(4:5),c(3),c(2),c(1)); 
disp('-----------------------------------------------------------------------------');

%-----------------------------------------------------------------------------------------------------------------------
%% LOADING DATA FROM MATFILE
%-----------------------------------------------------------------------------------------------------------------------

disp('Loading in raw data');

% Loading data
rec             = xlsread('Raw Files/USREC','USREC','b1549:b1956');
pmi             = xlsread('Raw Files/NAPM','NAPM','b376:b783');
cc              = xlsread('Raw Files/UMCSENT','UMCSENT','b16:b423');
ret             = xlsread('Raw Files/Data_bench','Data','k300:k707');
fed             = xlsread('Raw Files/Data_bench','Data','e300:e707');
tms             = xlsread('Raw Files/Data_bench','Data','d300:d707');
ip_growth       = xlsread('Raw Files/INDPRO','INDPRO','c756:c1157').* 100;
rawmacro        = xlsread('Raw Files/Panel2012','Macro','b165:ca573'); 
rawfinance      = xlsread('Raw Files/Panel2012','Financial','b167:cu575');
tcodemacro      = xlsread('Raw Files/Panel2012','Macro','b2:ca2');
tcodefinance    = xlsread('Raw Files/Panel2012','Financial','b3:cu3');

%-----------------------------------------------------------------------------------------------------------------------
%% ESTIMATING LATENT COMMON FACTORS BASED ON MACROECONOMIC AND FINANCIAL PANEL
%-----------------------------------------------------------------------------------------------------------------------

disp('Estimating latent common factors based on macroeconomic and financial variables');

% Combining macroeconomic and financial data
rawdata         = [rawmacro rawfinance];
tcode           = [tcodemacro tcodefinance];

% Setting up preliminaries
[nObs,nRows]    = size(rawdata);
mPanel          = nan(nObs,nRows);

% Looping over variables to transform
for iRow = 1:nRows

    if tcode(iRow) == 0 

        tcode(iRow) = 1; 

    end 

    % Transform variables according to transformation codes
    temp              = transx(rawdata(:,iRow),tcode(iRow));
    mPanel(:,iRow)    = temp;

end

% All data available from 1978M1 after transformation
mPanel = mPanel(2:end,:); 

% Estimating latent factors
[fhat,lambda,ehat] = pc_T(mPanel,20,'IC');

%-----------------------------------------------------------------------------------------------------------------------
%% GENERATING BENCHMARK VARIABLES
%-----------------------------------------------------------------------------------------------------------------------

disp('Generating benchmark variables');

% Choosing benchmark variables from the panel
ret2     = mPanel(:,98).*100;

fed2     = mPanel(:,83);
tms2     = mPanel(:,97);

%-----------------------------------------------------------------------------------------------------------------------
%% GENERATING FULL SET OF VARIABLES AND LAGS TO CHOOSE FROM
%-----------------------------------------------------------------------------------------------------------------------

disp('Generating full set of variables and lags to choose from');

% Preallocations
Lrec    = zeros(402,6);
Lpmi    = zeros(402,6);
Lcc     = zeros(402,6);
Ltms    = zeros(402,6);
Lfed    = zeros(402,6);
Lret    = zeros(402,6);
Lfhat1  = zeros(402,6);
Lfhat2  = zeros(402,6);
Lfhat3  = zeros(402,6);
Lfhat4  = zeros(402,6);
Lfhat5  = zeros(402,6);
Lfhat6  = zeros(402,6);
Lfhat7  = zeros(402,6);
Lfhat8  = zeros(402,6);
Lfhat9  = zeros(402,6);
Lfhat10 = zeros(402,6);
Lfhat11 = zeros(402,6);
Lfhat12 = zeros(402,6);
Lfhat13 = zeros(402,6);
Lfhat14 = zeros(402,6);
Lfhat15 = zeros(402,6);

% Looping over lag implementations
for h = 1:6

    Lrec(:,h)     = rec(7-h:end-h,1);
    Lpmi(:,h)     = pmi(7-h:end-h,1);
    Lcc(:,h)      = cc(7-h:end-h,1);
    Ltms(:,h)     = tms(7-h:end-h,1);
    Lfed(:,h)     = fed(7-h:end-h,1);
    Lret(:,h)     = ret(7-h:end-h,1);
    Lfhat1(:,h)   = fhat(7-h:end-h,1);
    Lfhat2(:,h)   = fhat(7-h:end-h,2);
    Lfhat3(:,h)   = fhat(7-h:end-h,3);
    Lfhat4(:,h)   = fhat(7-h:end-h,4);
    Lfhat5(:,h)   = fhat(7-h:end-h,5);
    Lfhat6(:,h)   = fhat(7-h:end-h,6);
    Lfhat7(:,h)   = fhat(7-h:end-h,7);
    Lfhat8(:,h)   = fhat(7-h:end-h,8);
    Lfhat9(:,h)   = fhat(7-h:end-h,9);
    Lfhat10(:,h)  = fhat(7-h:end-h,10);
    Lfhat11(:,h)  = fhat(7-h:end-h,11);
    Lfhat12(:,h)  = fhat(7-h:end-h,12);
    Lfhat13(:,h)  = fhat(7-h:end-h,13);
    Lfhat14(:,h)  = fhat(7-h:end-h,14);
    Lfhat15(:,h)  = fhat(7-h:end-h,15);

end

% Collecting lags of latent factors
Lfhat = [
    Lfhat1 ...
    Lfhat2 ...
    Lfhat3 ...
    Lfhat4 ...
    Lfhat5 ...
    Lfhat6 ...
    Lfhat7 ...
    Lfhat8 ...
    Lfhat9 ...
    Lfhat10 ...
    Lfhat11 ...
    Lfhat12 ...
    Lfhat13 ...
    Lfhat14 ...
    Lfhat15
];

% Recession one-period ahead (target variable)
frec = rec(7:end);

%-----------------------------------------------------------------------------------------------------------------------
%% SAVING RELEVANT VARIABLES
%-----------------------------------------------------------------------------------------------------------------------

save('matfiles/data.mat','frec','rec','pmi','cc','Lrec','Lpmi','Lcc','Ltms','Lfed','Lret','Lfhat',...
    'Lfhat1','Lfhat2','Lfhat3','Lfhat4','Lfhat5','Lfhat6','Lfhat7','Lfhat8','Lfhat9','Lfhat10','Lfhat11','Lfhat12',...
    'Lfhat13','Lfhat14','Lfhat15','tms','fed','ret','fhat','rawdata','tcode','mPanel','ip_growth');

%-----------------------------------------------------------------------------------------------------------------------
%% COMPUTING CODE RUN TIME
%-----------------------------------------------------------------------------------------------------------------------

tEnd = toc(tStart); rmpath('../');
fprintf('Runtime: %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
disp('Routine Completed');

%-----------------------------------------------------------------------------------------------------------------------
% END OF SCRIPT
%-----------------------------------------------------------------------------------------------------------------------