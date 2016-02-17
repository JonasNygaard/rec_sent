function [Lfhat1t, Lfhat2t, Lfhat3t, Lfhat4t, Lfhat5t, Lfhat6t, Lfhat7t, Lfhat8t, Lfhat9t, Lfhat10t,...
 Lfhat11t, Lfhat12t, Lfhat13t, Lfhat14t, Lfhat15t] = recursive_factor_construction(rawdata,tcode,nfac)

% MakeFactorsRecursively.m
%----------------------------------------------------------------------------------------------
% This function computes factors depending on the input data given. It is used for 
% the out-of-sample forecasting exercise. 
%
% Last modified: March 10, 2014
%----------------------------------------------------------------------------------------------

[T, N] = size(rawdata);
panel = nan(T,N);
for i=1:N;
      if tcode(i)==0; tcode(i)=1; end; 
       dum=transx(rawdata(:,i),tcode(i));
       panel(:,i) = dum;
end;
panel=panel(2:end,:); % All data available from 1978M1 after transformation
fhat = pc_T(panel,nfac,'Full');

nobs = size(fhat,1)-6;

Lfhat1t = zeros(nobs,6);
Lfhat2t = zeros(nobs,6);
Lfhat3t = zeros(nobs,6);
Lfhat4t = zeros(nobs,6);
Lfhat5t = zeros(nobs,6);
Lfhat6t = zeros(nobs,6);
Lfhat7t = zeros(nobs,6);
Lfhat8t = zeros(nobs,6);
Lfhat9t = zeros(nobs,6);
Lfhat10t = zeros(nobs,6);
Lfhat11t = zeros(nobs,6);
Lfhat12t = zeros(nobs,6);
Lfhat13t = zeros(nobs,6);
Lfhat14t = zeros(nobs,6);
Lfhat15t = zeros(nobs,6);
for h = 1:6
      Lfhat1t(:,h) = fhat(7-h:end-h,1);
      Lfhat2t(:,h) = fhat(7-h:end-h,2);
      Lfhat3t(:,h) = fhat(7-h:end-h,3);
      Lfhat4t(:,h) = fhat(7-h:end-h,4);
      Lfhat5t(:,h) = fhat(7-h:end-h,5);
      Lfhat6t(:,h) = fhat(7-h:end-h,6);
      Lfhat7t(:,h) = fhat(7-h:end-h,7);
      Lfhat8t(:,h) = fhat(7-h:end-h,8);
      Lfhat9t(:,h) = fhat(7-h:end-h,9);
      Lfhat10t(:,h) = fhat(7-h:end-h,10);
      Lfhat11t(:,h) = fhat(7-h:end-h,11);
      Lfhat12t(:,h) = fhat(7-h:end-h,12);
      Lfhat13t(:,h) = fhat(7-h:end-h,13);
      Lfhat14t(:,h) = fhat(7-h:end-h,14);
      Lfhat15t(:,h) = fhat(7-h:end-h,15);
end
end

%----------------------------------------------------------------------------------------------
% END OF FUNCTION
%----------------------------------------------------------------------------------------------