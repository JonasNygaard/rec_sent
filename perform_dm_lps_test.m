function [DMvalue,pval] = DMtest_lps(lps1,lps2)

% DMtest.m
%----------------------------------------------------------------------------------------------
% This function computes the Diebold-Mariano test statistics for two
% competing forecast from non-nested models. 
%
% Last modified: March 19, 2014
%----------------------------------------------------------------------------------------------

diff_lps    = lps1 - lps2;
iota        = ones(size(lps1,1),1);
res         = nwest(diff_lps,iota,3);

% [DMvalue,sbv] = olsnwhh(diff_lps,iota,3,'NW',0);
% res_tms_stat        = nwest(ip_growth,[ones(402,1) Ltms(:,4)],3);

DMvalue     = res.beta;
pval        = 1-normcdf(res.tstat,0,1);

end

%----------------------------------------------------------------------------------------------
% END OF FUNCTION
%----------------------------------------------------------------------------------------------