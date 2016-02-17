function [testval,pval] = perform_auc_test(actual,auc1,auc2,forecast1,forecast2)

% dAUCtest
%----------------------------------------------------------------------------------------------
% This function computes the test for a difference in two AUC values. 
%
% Last modified: March 19, 2014
%----------------------------------------------------------------------------------------------

indx = find(actual == 1);
n1 = size(actual(indx,1),1);
r1 = corr(forecast1(indx),forecast2(indx));
indx = find(actual == 0);
r0 = corr(forecast1(indx),forecast2(indx));
n0 = size(actual(indx,1),1);
r = (r1+r0)/2;

Q1 = auc1 / (2-auc1);
Q2 = (2*auc1^2) / (1+auc1);
se1 = sqrt((1/(n0*n1))*(auc1*(1-auc1)+(n1-1)*(Q1-auc1^2)+(n0-1)*(Q2-auc1^2)));

Q1 = auc2 / (2-auc2);
Q2 = (2*auc2^2) / (1+auc2);
se2 = sqrt((1/(n0*n1))*(auc2*(1-auc2)+(n1-1)*(Q1-auc2^2)+(n0-1)*(Q2-auc2^2)));

testval = (auc1 - auc2) / sqrt(se1^2 + se2^2 - 2*r*se1*se2);
pval = 1-normcdf(testval,0,1);

end

%----------------------------------------------------------------------------------------------
% END OF FUNCTION
%----------------------------------------------------------------------------------------------