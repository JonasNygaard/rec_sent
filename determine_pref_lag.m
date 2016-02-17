%% determine_pref_lags.m
%-----------------------------------------------------------------------------------------------------------------------
%
%   This script determines the preferred lag length for each of the variables that we consider in the paper: 
%   "Forecasting US Recessions: The Role of Sentiment".  
%
%   --------------------------------
%   Last modified: December, 2015
%   --------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------

disp('Running the determine_pref_lag.m script');

%-----------------------------------------------------------------------------------------------------------------------
%% SETTING UP PRELIMINARIES
%-----------------------------------------------------------------------------------------------------------------------

% Preallocations
prefdata    = [];
strall      = {};
strm        = {};

%-----------------------------------------------------------------------------------------------------------------------
%% DETERMINING PREFERRED LAG LENGTHS
%-----------------------------------------------------------------------------------------------------------------------

% Finding preferred lag length for PMI
str = {'PMI$_{t-1}$','PMI$_{t-2}$','PMI$_{t-3}$','PMI$_{t-4}$','PMI$_{t-5}$','PMI$_{t-6}$'};
bic = zeros(6,1);
for h = 1:6

    res         = statProbit(frec,Lpmi(:,h));
    bic(h,1)    = res.bic;

end
[~,min_ind] = min(bic);
strall{1}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lpmi(:,min_ind)];

% Finding preferred lag length for CC
str = {'CC$_{t-1}$','CC$_{t-2}$','CC$_{t-3}$','CC$_{t-4}$','CC$_{t-5}$','CC$_{t-6}$'};
bic = zeros(6,1);
for h = 1:6

    res         = statProbit(frec,Lcc(:,h));
    bic(h,1)    = res.bic;

end
[~,min_ind] = min(bic);
strall{2}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lcc(:,min_ind)];

% Finding preferred lag length for TS
str = {'TS$_{t-1}$','TS$_{t-2}$','TS$_{t-3}$','TS$_{t-4}$','TS$_{t-5}$','TS$_{t-6}$'};
bic = zeros(6,1);
for h = 1:6

    res = statProbit(frec,Ltms(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{3}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Ltms(:,min_ind)];

% Finding preferred lag length for FFR
str = {'FFR$_{t-1}$','FFR$_{t-2}$','FFR$_{t-3}$','FFR$_{t-4}$','FFR$_{t-5}$','FFR$_{t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfed(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{4}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfed(:,min_ind)];

% Finding preferred lag length for RET
str = {'RET$_{t-1}$','RET$_{t-2}$','RET$_{t-3}$','RET$_{t-4}$','RET$_{t-5}$','RET$_{t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lret(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{5}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lret(:,min_ind)];

% Finding preferred lag length for F1
str = {'$\hat{f}_{1,t-1}$','$\hat{f}_{1,t-2}$','$\hat{f}_{1,t-3}$','$\hat{f}_{1,t-4}$',...
'$\hat{f}_{1,t-5}$','$\hat{f}_{1,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat1(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{6}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat1(:,min_ind)];

% Finding preferred lag length for F2
str = {'$\hat{f}_{2,t-1}$','$\hat{f}_{2,t-2}$','$\hat{f}_{2,t-3}$','$\hat{f}_{2,t-4}$',...
'$\hat{f}_{2,t-5}$','$\hat{f}_{2,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat2(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{7}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat2(:,min_ind)];

% Finding preferred lag length for F3
str = {'$\hat{f}_{3,t-1}$','$\hat{f}_{3,t-2}$','$\hat{f}_{3,t-3}$','$\hat{f}_{3,t-4}$',...
'$\hat{f}_{3,t-5}$','$\hat{f}_{3,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat3(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{8}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat3(:,min_ind)];

% Finding preferred lag length for F4
str = {'$\hat{f}_{4,t-1}$','$\hat{f}_{4,t-2}$','$\hat{f}_{4,t-3}$','$\hat{f}_{4,t-4}$',...
'$\hat{f}_{4,t-5}$','$\hat{f}_{4,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat4(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{9}   = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat4(:,min_ind)];

% Finding preferred lag length for F5
str = {'$\hat{f}_{5,t-1}$','$\hat{f}_{5,t-2}$','$\hat{f}_{5,t-3}$','$\hat{f}_{5,t-4}$',...
'$\hat{f}_{5,t-5}$','$\hat{f}_{5,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat5(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{10}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat5(:,min_ind)];

% Finding preferred lag length for F6
str = {'$\hat{f}_{6,t-1}$','$\hat{f}_{6,t-2}$','$\hat{f}_{6,t-3}$','$\hat{f}_{6,t-4}$',...
'$\hat{f}_{6,t-5}$','$\hat{f}_{6,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat6(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{11}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat6(:,min_ind)];

% Finding preferred lag length for F7
str = {'$\hat{f}_{7,t-1}$','$\hat{f}_{7,t-2}$','$\hat{f}_{7,t-3}$','$\hat{f}_{7,t-4}$',...
'$\hat{f}_{7,t-5}$','$\hat{f}_{7,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat7(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{12}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat7(:,min_ind)];

% Finding preferred lag length for F8
str = {'$\hat{f}_{8,t-1}$','$\hat{f}_{8,t-2}$','$\hat{f}_{8,t-3}$','$\hat{f}_{8,t-4}$',...
'$\hat{f}_{8,t-5}$','$\hat{f}_{8,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat8(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{13}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat8(:,min_ind)];

% Finding preferred lag length for F9
str = {'$\hat{f}_{9,t-1}$','$\hat{f}_{9,t-2}$','$\hat{f}_{9,t-3}$','$\hat{f}_{9,t-4}$',...
'$\hat{f}_{9,t-5}$','$\hat{f}_{9,t-6}$'};
bic = zeros(6,1);
for h = 1:6

    res = statProbit(frec,Lfhat9(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{14}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat9(:,min_ind)];

% Finding preferred lag length for F10
str = {'$\hat{f}_{10,t-1}$','$\hat{f}_{10,t-2}$','$\hat{f}_{10,t-3}$','$\hat{f}_{10,t-4}$',...
'$\hat{f}_{10,t-5}$','$\hat{f}_{10,t-6}$'};
bic = zeros(6,1);
for h = 1:6

    res = statProbit(frec,Lfhat10(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{15}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat10(:,min_ind)];

% Finding preferred lag length for F11
str = {'$\hat{f}_{11,t-1}$','$\hat{f}_{11,t-2}$','$\hat{f}_{11,t-3}$','$\hat{f}_{11,t-4}$',...
'$\hat{f}_{11,t-5}$','$\hat{f}_{11,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat11(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{16}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat11(:,min_ind)];

% Finding preferred lag length for F12
str = {'$\hat{f}_{12,t-1}$','$\hat{f}_{12,t-2}$','$\hat{f}_{12,t-3}$','$\hat{f}_{12,t-4}$',...
'$\hat{f}_{12,t-5}$','$\hat{f}_{12,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat12(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{17}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat12(:,min_ind)];

% Finding preferred lag length for F13
str = {'$\hat{f}_{13,t-1}$','$\hat{f}_{13,t-2}$','$\hat{f}_{13,t-3}$','$\hat{f}_{13,t-4}$',...
'$\hat{f}_{13,t-5}$','$\hat{f}_{13,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat13(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{18}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat13(:,min_ind)];

% Finding preferred lag length for F14
str = {'$\hat{f}_{14,t-1}$','$\hat{f}_{14,t-2}$','$\hat{f}_{14,t-3}$','$\hat{f}_{14,t-4}$',...
'$\hat{f}_{14,t-5}$','$\hat{f}_{14,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat14(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{19}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat14(:,min_ind)];

% Finding preferred lag length for F15
str = {'$\hat{f}_{15,t-1}$','$\hat{f}_{15,t-2}$','$\hat{f}_{15,t-3}$','$\hat{f}_{15,t-4}$',...
'$\hat{f}_{15,t-5}$','$\hat{f}_{15,t-6}$'};
bic = zeros(6,1);
for h = 1:6
    
    res = statProbit(frec,Lfhat15(:,h));
    bic(h,1) = res.bic;

end
[~,min_ind] = min(bic);
strall{20}  = str{min_ind};
strm        = [strm str];
prefdata    = [prefdata Lfhat15(:,min_ind)];

% Set string vector for tables
strm = [strm {'$\pi_{t-1}$','y$_{t-3}$'}];

%-----------------------------------------------------------------------------------------------------------------------
%% END OF FUNCTION
%-----------------------------------------------------------------------------------------------------------------------