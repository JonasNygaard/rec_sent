%% main_cem.m
%-----------------------------------------------------------------------------------------------------------------------
%
%   This scripts contain the main code for the paper "Forecasting US Recessions: The Role of Sentiment". 
%
%   The script imports data from the Data/Output/ folder, which have been prepared using the data making scripts stored
%   in the Data/ folder. The Data/ folder contains a Markdown file with further notes on the data. Using these data, 
%   this script compute all results presented in the paper following the procedures described in more detail in the 
%   paper. All results are stored in LaTeX tables and encapsulated post script (eps) figures in the TeX/ folder and 
%   used directly as input to the main TeX file for the paper. The script calls all auxiliary functions and scripts 
%   from the root directory when needed along the way. 
%    
%   -------------------------------
%   Last modified: December, 2015
%   -------------------------------
%
%-----------------------------------------------------------------------------------------------------------------------

clear; clc; close all; tStart = tic; format shortg; c = clock; 
disp('-----------------------------------------------------------------------------');
disp('Main script for the paper "Forecasting US Recessions: The Role of Sentiment" ');
fprintf('Code initiated at %.0f:%.0f on %.0f / %0.f - %0.f \n',c(4:5),c(3),c(2),c(1)); 
disp('-----------------------------------------------------------------------------');

%-----------------------------------------------------------------------------------------------------------------------
%% LOADING AND PREPARING DATA
%-----------------------------------------------------------------------------------------------------------------------
%{
    Data are constructed in the Data/ folder using the Matlab files stored therein. Here we load the processed data 
    from mat files and set it up for our empirical analysis and sample period. 
%}
%-----------------------------------------------------------------------------------------------------------------------

disp('Loading and preparing data');

% Loading data files
load('Data/matfiles/data.mat');

% Run script to determine preferred lag for each predictor variable
run('determine_pref_lag.m');

% Creating datenum vector
datenum_indx    = datenum({'01-Jan-1978';'02-Jan-2012'});
datevec_indx    = datevec(datenum_indx(1):1:datenum_indx(2));
datevec_indx    = datevec_indx(datevec_indx(2:end,2) ~= datevec_indx(1:end-1,2),1:3);
datenum_indx    = datenum(datevec_indx);

%-----------------------------------------------------------------------------------------------------------------------
%% PLOTTING SENTIMENT INDEXES AGAINS NBER RECESSION DATES
%-----------------------------------------------------------------------------------------------------------------------

disp('Plotting sentiment indexes against NBER recession dates');

% Collecting sentiment variables
sent    = [pmi cc];

% Plotting sentiment indexes against NBER recession periods
figure;
hold on 
b1 = bar(datenum_indx,rec.*119.5);
b1.EdgeColor    = [0.8 0.8 0.8];
b1.FaceColor    = [0.8 0.8 0.8];
b1.BarWidth     = 1;
b1.ShowBaseLine = 'off';
set(get(get(b1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
p1 = plot(datenum_indx,sent);
p1(1).LineStyle = '-';  p1(1).Color = 'k';
p1(2).LineStyle = '-.'; p1(2).Color = 'k';
p1(1).LineWidth = 1.1;  p1(2).LineWidth = 1.1;
hold off
axis([-inf inf 0 120]);
datetick('x','yyyy');
leg = legend('PMI$_{t}$','CC$_{t}$');
set(leg,'Position',[0.7  0.2 0.1 0.03],'Box','off','Interpreter','latex');
box on
set(gcf,'PaperUnits','Centimeters','PaperSize',[23 9],'PaperPosition',[0.5 0.5 22 8]);
print(gcf,'-depsc','TeX/Figures/Figure1.eps');

%-----------------------------------------------------------------------------------------------------------------------
%% COMPUTING SUMMARY STATISTICS FOR PREDICTOR VARIBLES
%-----------------------------------------------------------------------------------------------------------------------

disp('Computing summary statistics for predictor variables');

% Collecting predictor variables
sumdata = [pmi cc tms fed ret fhat];

% Computing correlation coefficients
cormat  = corr(sumdata);

% Creating LaTex Table for correlation coefficient
fid = fopen('TeX/Tables/cem_correlations.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s \\\\\\midrule\n','Variable',strm{1},strm{7},'Variable',strm{1},strm{7});
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{13},cormat(3,1),cormat(3,2),strm{67},cormat(12,1),cormat(12,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{19},cormat(4,1),cormat(4,2),strm{73},cormat(13,1),cormat(13,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{25},cormat(5,1),cormat(5,2),strm{79},cormat(14,1),cormat(14,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{31},cormat(6,1),cormat(6,2),strm{85},cormat(15,1),cormat(15,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{37},cormat(7,1),cormat(7,2),strm{91},cormat(16,1),cormat(16,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{43},cormat(8,1),cormat(8,2),strm{97},cormat(17,1),cormat(17,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{49},cormat(9,1),cormat(9,2),strm{103},cormat(18,1),cormat(18,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{55},cormat(10,1),cormat(10,2),strm{109},cormat(19,1),cormat(19,2));
fprintf(fid,'%s & %.2f & %.2f & %s & %.2f & %.2f \\\\\n',strm{61},cormat(11,1),cormat(11,2),strm{115},cormat(20,1),cormat(20,2));
fprintf(fid,'\n');
fclose(fid);

%-----------------------------------------------------------------------------------------------------------------------
%% ESTIMATING SINGLE-FACTOR PROBIT MODELS USING PREFERRED LAGS FOR EACH PREDICTOR
%-----------------------------------------------------------------------------------------------------------------------

disp('Estimating single-factor Probit models using preferred lags for each predictor');

% Preallocations
nVar        = size(prefdata,2);
nObs        = size(frec,1);
pseudo_r2   = NaN(nVar,1);
bic         = NaN(nVar,1);
auc         = NaN(nVar,1);
aucp        = NaN(nVar,1);
likhood     = NaN(nVar,1);
lps         = NaN(nVar,1);
prob        = NaN(nObs,nVar);

% Estimate models one at a time
for iVar = 1:nVar

    res                 = statProbit(frec,prefdata(:,iVar));
    pseudo_r2(iVar,1)   = res.r2est;
    bic(iVar,1)         = res.bic;
    auc(iVar,1)         = res.auc;
    aucp(iVar,1)        = res.aucp;
    likhood(iVar,1)     = res.likhood;
    lps(iVar,1)         = res.lps;
    prob(:,iVar)        = res.yhat;

end

% Sort variables according to BIC criterion
[~,bic_indx]    = sort(bic(:,1));

% Creating LaTeX Table
fid = fopen('TeX/Tables/cem_univariate_models.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n','Rank','Variable','pseudo-R$^{2}$','BIC','LPS','AUC','p-val');
for k=1:size(prefdata,2)
    fprintf(fid,'%.0f & %s & %.2f & %.2f & %.2f & %.2f & (%.2f) \\\\\n',k,strall{bic_indx(k)},...   
    pseudo_r2(bic_indx(k)),bic(bic_indx(k)),lps(bic_indx(k)),auc(bic_indx(k)),aucp(bic_indx(k)));
end
fprintf(fid,'\n');
fclose(fid);

%-----------------------------------------------------------------------------------------------------------------------
%% ESTIMATING VARIOUS IN-SAMPLE PROBIT MODELS
%-----------------------------------------------------------------------------------------------------------------------

disp('Estimating various in-sample probit models');

% Setting the publication lag
Llrec               = Lrec(:,3);

% Classic predictor models
res_tms_stat        = statProbit(frec,Ltms(:,6));
res_class_stat      = statProbit(frec,[Ltms(:,6) Lfed(:,6) Lret(:,[2 4 6])]);
res_class_dyn       = statProbit(frec,[Ltms(:,3) Lret(:,[2 4]) Llrec]);
res_class_ar        = dynProbit(frec,[Ltms(:,6) Lret(:,2)]);
res_class_dynar     = dynProbit(frec,[Ltms(:,6) Lret(:,2) Llrec]);

% Sentiment predictor models
res_pmi_stat        = statProbit(frec,Lpmi(:,1));
res_cc_stat         = statProbit(frec,Lcc(:,1));
res_sent_stat       = statProbit(frec,[Lpmi(:,1) Lcc(:,[1 6])]);
res_sent_dyn        = statProbit(frec,[Lpmi(:,1) Lcc(:,[1 5]) Llrec]);
res_sent_ar         = dynProbit(frec,[Lpmi(:,[1 2 4]) Lcc(:,[1 6])]);
res_sent_dynar      = dynProbit(frec,[Lpmi(:,[1 2 4]) Lcc(:,[1 6]) Llrec]);

% Sentiment and classic predictor models
res_classent_stat   = statProbit(frec,[Lpmi(:,1) Lcc(:,1) Ltms(:,6) Lret(:,[2 4])]);
res_classent_dyn    = statProbit(frec,[Lpmi(:,1) Lcc(:,1) Ltms(:,6) Lret(:,[2 4]) Llrec]);
res_classent_ar     = dynProbit(frec,[Lpmi(:,6) Lcc(:,[1 5]) Ltms(:,6) Lret(:,2)]);
res_classent_dynar  = dynProbit(frec,[Lpmi(:,5) Lcc(:,[1 5]) Ltms(:,6) Lret(:,2) Llrec]);

% Sentiment and factor predictor models
res_facsent_stat    = statProbit(frec,[Lpmi(:,1) Lcc(:,1) Lfhat1(:,2) Lfhat6(:,6) Lfhat14(:,1)]);
res_facsent_dyn     = statProbit(frec,[Lpmi(:,1) Lcc(:,1) Lfhat1(:,2) Lfhat6(:,6) Lfhat14(:,1) Llrec]);
res_facsent_ar      = dynProbit(frec,[Lpmi(:,1) Lfhat1(:,2) Lfhat2(:,2) Lfhat6(:,6) Lfhat14(:,1)]);
res_facsent_dynar   = dynProbit(frec,[Lpmi(:,3) Lfhat1(:,2) Lfhat2(:,2) Lfhat6(:,6) Lfhat14(:,1) Llrec]);

% Creating LaTeX table for static probit models
fid = fopen('TeX/Tables/cem_in_sample_results.tex','w');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %s \\\\\\midrule\n','Variable','Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7');
fprintf(fid,'%s & %s & %s & %.2f & %s & %.2f & %.2f & %.2f \\\\\n',strm{1},'','',res_pmi_stat.beta(2),'',res_sent_stat.beta(2),res_classent_stat.beta(2),res_facsent_stat.beta(2));
fprintf(fid,'%s & %s & %s & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) \\\\\n','','','',res_pmi_stat.stdb(2),'',res_sent_stat.stdb(2),res_classent_stat.stdb(2),res_facsent_stat.stdb(2));
fprintf(fid,'%s & %s & %s & %s & %.2f & %.2f & %.2f & %.2f \\\\\n',strm{7},'','','',res_cc_stat.beta(2),res_sent_stat.beta(3),res_classent_stat.beta(3),res_facsent_stat.beta(3));
fprintf(fid,'%s & %s & %s & %s & (%.2f) & (%.2f) & (%.2f) & (%.2f) \\\\\n','','','','',res_cc_stat.stdb(2),res_sent_stat.stdb(3),res_classent_stat.stdb(3),res_facsent_stat.stdb(3));
fprintf(fid,'%s & %s & %s & %s & %s & %.2f & %s & %s \\\\\n',strm{12},'','','','',res_sent_stat.beta(4),'','');
fprintf(fid,'%s & %s & %s & %s & %s & (%.2f) & %s & %s \\\\\n','','','','','',res_sent_stat.stdb(4),'','');
fprintf(fid,'%s & %.2f & %.2f & %s & %s & %s & %.2f & %s \\\\\n',strm{18},res_tms_stat.beta(2),res_class_stat.beta(2),'','','',res_classent_stat.beta(4),'');
fprintf(fid,'%s & (%.2f) & (%.2f) & %s & %s & %s & (%.2f) & %s \\\\\n','',res_tms_stat.stdb(2),res_class_stat.stdb(2),'','','',res_classent_stat.stdb(4),'');
fprintf(fid,'%s & %s & %.2f & %s & %s & %s & %s & %s \\\\\n',strm{24},'',res_class_stat.beta(3),'','','','','');
fprintf(fid,'%s & %s & (%.2f) & %s & %s & %s & %s & %s \\\\\n','','',res_class_stat.stdb(3),'','','','','');
fprintf(fid,'%s & %s & %.2f & %s & %s & %s & %.2f & %s \\\\\n',strm{26},'',res_class_stat.beta(4),'','','',res_classent_stat.beta(5),'');
fprintf(fid,'%s & %s & (%.2f) & %s & %s & %s & (%.2f) & %s \\\\\n','','',res_class_stat.stdb(4),'','','',res_classent_stat.stdb(5),'');
fprintf(fid,'%s & %s & %.2f & %s & %s & %s & %.2f & %s \\\\\n',strm{28},'',res_class_stat.beta(5),'','','',res_classent_stat.beta(6),'');
fprintf(fid,'%s & %s & (%.2f) & %s & %s & %s & (%.2f) & %s \\\\\n','','',res_class_stat.stdb(5),'','','',res_classent_stat.stdb(6),'');
fprintf(fid,'%s & %s & %.2f & %s & %s & %s & %s & %s \\\\\n',strm{30},'',res_class_stat.beta(6),'','','','','');
fprintf(fid,'%s & %s & (%.2f) & %s & %s & %s & %s & %s \\\\\n','','',res_class_stat.stdb(6),'','','','','');
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %.2f \\\\\n',strm{32},'','','','','','',res_facsent_stat.beta(4));
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & (%.2f) \\\\\n','','','','','','','',res_facsent_stat.stdb(4));
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %.2f \\\\\n',strm{66},'','','','','','',res_facsent_stat.beta(5));
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & (%.2f) \\\\\n','','','','','','','',res_facsent_stat.stdb(5));
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & %.2f \\\\\n',strm{112},'','','','','','',res_facsent_stat.beta(6));
fprintf(fid,'%s & %s & %s & %s & %s & %s & %s & (%.2f) \\\\\n','','','','','','','',res_facsent_stat.stdb(6));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','log-$\mathcal{L}$',...
res_tms_stat.likhood,res_class_stat.likhood,res_pmi_stat.likhood,res_cc_stat.likhood,res_sent_stat.likhood,res_classent_stat.likhood,res_facsent_stat.likhood);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','pseudo-R$^{2}$',...
res_tms_stat.r2est,res_class_stat.r2est,res_pmi_stat.r2est,res_cc_stat.r2est,res_sent_stat.r2est,res_classent_stat.r2est,res_facsent_stat.r2est);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','BIC',...
res_tms_stat.bic,res_class_stat.bic,res_pmi_stat.bic,res_cc_stat.bic,res_sent_stat.bic,res_classent_stat.bic,res_facsent_stat.bic);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','LPS',...
res_tms_stat.lps,res_class_stat.lps,res_pmi_stat.lps,res_cc_stat.lps,res_sent_stat.lps,res_classent_stat.lps,res_facsent_stat.lps);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','CR50\%',...
res_tms_stat.cr50,res_class_stat.cr50,res_pmi_stat.cr50,res_cc_stat.cr50,res_sent_stat.cr50,res_classent_stat.cr50,res_facsent_stat.cr50);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','CE50\%',...
res_tms_stat.ce50,res_class_stat.ce50,res_pmi_stat.ce50,res_cc_stat.ce50,res_sent_stat.ce50,res_classent_stat.ce50,res_facsent_stat.ce50);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','AUC',...
res_tms_stat.auc,res_class_stat.auc,res_pmi_stat.auc,res_cc_stat.auc,res_sent_stat.auc,res_classent_stat.auc,res_facsent_stat.auc);
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) & (%.2f) \\\\\n','p-val',...
res_tms_stat.aucp,res_class_stat.aucp,res_pmi_stat.aucp,res_cc_stat.aucp,res_sent_stat.aucp,res_classent_stat.aucp,res_facsent_stat.aucp);
fprintf(fid,'\n');
fclose(fid);

% Creating LaTeX Table for robustness check on dynamic models
fid = fopen('TeX/Tables/cem_dynamic_probit_models.tex','w');
fprintf(fid,'%s & %s & %s & %s \\\\\\midrule\n','Variable','Dynamic 1','Dynamic 2','Dynamic 3');
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n',strm{1},res_sent_dyn.beta(2),res_classent_dyn.beta(2),res_facsent_dyn.beta(2));
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) \\\\\n','',res_sent_dyn.stdb(2),res_classent_dyn.stdb(2),res_facsent_dyn.stdb(2));
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n',strm{7},res_sent_dyn.beta(3),res_classent_dyn.beta(3),res_facsent_dyn.beta(3));
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) \\\\\n','',res_sent_dyn.stdb(3),res_classent_dyn.stdb(3),res_facsent_dyn.stdb(3));
fprintf(fid,'%s & %.2f & %s & %s \\\\\n',strm{11},res_sent_dyn.beta(4),'','');
fprintf(fid,'%s & (%.2f) & %s & %s \\\\\n','',res_sent_dyn.stdb(4),'','');
fprintf(fid,'%s & %s & %.2f & %s \\\\\n',strm{18},'',res_classent_dyn.beta(4),'');
fprintf(fid,'%s & %s & (%.2f) & %s \\\\\n','','',res_classent_dyn.stdb(4),'');
fprintf(fid,'%s & %s & %.2f & %s \\\\\n',strm{26},'',res_classent_dyn.beta(5),'');
fprintf(fid,'%s & %s & (%.2f) & %s \\\\\n','','',res_classent_dyn.stdb(5),'');
fprintf(fid,'%s & %s & %.2f & %s \\\\\n',strm{28},'',res_classent_dyn.beta(6),'');
fprintf(fid,'%s & %s & (%.2f) & %s \\\\\n','','',res_classent_dyn.stdb(6),'');
fprintf(fid,'%s & %s & %s & %.2f \\\\\n',strm{32},'','',res_facsent_dyn.beta(4));
fprintf(fid,'%s & %s & %s & (%.2f) \\\\\n','','','',res_facsent_dyn.stdb(4));
fprintf(fid,'%s & %s & %s & %.2f \\\\\n',strm{66},'','',res_facsent_dyn.beta(5));
fprintf(fid,'%s & %s & %s & (%.2f) \\\\\n','','','',res_facsent_dyn.stdb(5));
fprintf(fid,'%s & %s & %s & %.2f \\\\\n',strm{112},'','',res_facsent_dyn.beta(6));
fprintf(fid,'%s & %s & %s & (%.2f) \\\\\n','','','',res_facsent_dyn.stdb(6));
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n',strm{122},res_sent_dyn.beta(5),res_classent_dyn.beta(7),res_facsent_dyn.beta(7));
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) \\\\\n','',res_sent_dyn.stdb(5),res_classent_dyn.stdb(7),res_facsent_dyn.stdb(7));
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n','log-$\mathcal{L}$',res_sent_dyn.likhood,res_classent_dyn.likhood,res_facsent_dyn.likhood);
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n','pseudo-R$^{2}$',res_sent_dyn.r2est,res_classent_dyn.r2est,res_facsent_dyn.r2est);
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n','BIC',res_sent_dyn.bic,res_classent_dyn.bic,res_facsent_dyn.bic);
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n','LPS',res_sent_dyn.lps,res_classent_dyn.lps,res_facsent_dyn.lps);
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n','CR50\%',res_sent_dyn.cr50,res_classent_dyn.cr50,res_facsent_dyn.cr50);
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n','CE50\%',res_sent_dyn.ce50,res_classent_dyn.ce50,res_facsent_dyn.ce50);
fprintf(fid,'%s & %.2f & %.2f & %.2f \\\\\n','AUC',res_sent_dyn.auc,res_classent_dyn.auc,res_facsent_dyn.auc);
fprintf(fid,'%s & (%.2f) & (%.2f) & (%.2f) \\\\\n','p-val',res_sent_dyn.aucp,res_classent_dyn.aucp,res_facsent_dyn.aucp);
fprintf(fid,'\n');
fclose(fid);

% Computing test for difference in AUC
[dauc_class_classentIS,daucp_class_classentIS]      = perform_auc_test(frec,res_classent_stat.auc,res_class_stat.auc,res_classent_stat.yhat,res_class_stat.yhat);
[dauc_class_sentIS,daucp_class_sentIS]              = perform_auc_test(frec,res_sent_stat.auc,res_class_stat.auc,res_sent_stat.yhat,res_class_stat.yhat);
[dauc_class_facsentIS,daucp_class_facsentIS]        = perform_auc_test(frec,res_facsent_stat.auc,res_class_stat.auc,res_facsent_stat.yhat,res_class_stat.yhat);
[dauc_sent_classentIS,daucp_sent_classentIS]        = perform_auc_test(frec,res_classent_stat.auc,res_sent_stat.auc,res_classent_stat.yhat,res_sent_stat.yhat);
[dauc_sent_facsentIS,daucp_sent_facsentIS]          = perform_auc_test(frec,res_facsent_stat.auc,res_sent_stat.auc,res_facsent_stat.yhat,res_sent_stat.yhat);
[dauc_classent_facsentIS,daucp_classent_facsentIS]  = perform_auc_test(frec,res_facsent_stat.auc,res_classent_stat.auc,res_facsent_stat.yhat,res_classent_stat.yhat);

%-----------------------------------------------------------------------------------------------------------------------
%% PLOTTING THE IN-SAMPLE FIT FOR SELECTED MODELS
%-----------------------------------------------------------------------------------------------------------------------

disp('PLotting the in-sample fit for selected models');

% Create text string
model_str = {
    'Model 2'
    'Model 5'
    'Model 6'
    'Model 7'
};

% Collect implied recession probabilities for selected models
prob_rec = [
    res_class_stat.yhat ...
    res_sent_stat.yhat ...
    res_classent_stat.yhat ...
    res_facsent_stat.yhat
];

% Plotting implied recession probabilities
figure;
for iFig = 1:4

    subplot(2,2,iFig);
    hold on

    b1 = bar(datenum_indx(7:end,:),frec);
    b1.EdgeColor = [0.8 0.8 0.8];
    b1.FaceColor = [0.8 0.8 0.8];
    b1.BarWidth  = 1;
    p1 = plot(datenum_indx(7:end,:),prob_rec(:,iFig));
    p1.LineStyle = '-'; p1.Color = 'k';
    thinl = line([datenum(1978,1,1) ; datenum(2012,1,1)],[0.5 ; 0.5],'LineStyle',':','Color','k');
    set(get(get(thinl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    datetick('x','yyyy');
    axis([-inf inf 0 1.004]);
    title(model_str{iFig},'FontSize',8);
    xlabel('Time');
    ylabel('Probability');
    box on

end
set(gcf,'PaperUnits','Centimeters','PaperSize',[21 15],'PaperPosition',[0.5 0.5 20 14]);
print(gcf,'-depsc','TeX/Figures/Figure2.eps');

%-----------------------------------------------------------------------------------------------------------------------
%% TESTING FOR PARAMETER INSTABILITY IN VARIOUS IN-SAMPLE PROBIT MODELS
%-----------------------------------------------------------------------------------------------------------------------

disp('Testing for parameter instability in various in-sample probit models');

% Setting trimming rule
trimr   = 0.25;

% Term-spread predictor models
[LM_tms_stat,cval_tms_stat]             = StabilityTestStat(frec,Ltms(:,6),trimr);
[LM_tms_dyn,cval_tms_dyn]               = StabilityTestStat(frec,[Ltms(:,6) Llrec],trimr);
[LM_tms_ar,cval_tms_ar]                 = StabilityTestDyn(frec,Ltms(:,6),trimr);
[LM_tms_dynar,cval_tms_dynar]           = StabilityTestDyn(frec,[Ltms(:,6) Llrec],trimr);

% PMI predictor models
[LM_pmi_stat,cval_pmi_stat]             = StabilityTestStat(frec,Lpmi(:,1),trimr);
[LM_pmi_dyn,cval_pmi_dyn]               = StabilityTestStat(frec,[Lpmi(:,1) Llrec],trimr);
[LM_pmi_ar,cval_pmi_ar]                 = StabilityTestDyn(frec,Lpmi(:,1),trimr);
[LM_pmi_dynar,cval_pmi_dynar]           = StabilityTestDyn(frec,[Lpmi(:,1) Llrec],trimr);

% CC predictor models
[LM_cc_stat,cval_cc_stat]               = StabilityTestStat(frec,Lcc(:,1),trimr);
[LM_cc_dyn,cval_cc_dyn]                 = StabilityTestStat(frec,[Lcc(:,1) Llrec],trimr);
[LM_cc_ar,cval_cc_ar]                   = StabilityTestDyn(frec,Lcc(:,1),trimr);
[LM_cc_dynar,cval_cc_dynar]             = StabilityTestDyn(frec,[Lcc(:,1) Llrec],trimr);

% Classic predictor models
[LM_class_stat,cval_class_stat]         = StabilityTestStat(frec,[Ltms(:,6) Lfed(:,6) Lret(:,[2 4 6])],trimr);
[LM_class_dyn,cval_class_dyn]           = StabilityTestStat(frec,[Ltms(:,3) Lret(:,[2 4]) Llrec],trimr);
[LM_class_ar,cval_class_ar]             = StabilityTestDyn(frec,[Ltms(:,6) Lret(:,2)],trimr);
[LM_class_dynar,cval_class_dynar]       = StabilityTestDyn(frec,[Ltms(:,6) Lret(:,2) Llrec],trimr);

% Sentiment predictor models
[LM_sent_stat,cval_sent_stat]           = StabilityTestStat(frec,[Lpmi(:,1) Lcc(:,[1 6])],trimr);
[LM_sent_dyn,cval_sent_dyn]             = StabilityTestStat(frec,[Lpmi(:,1) Lcc(:,[1 5]) Llrec],trimr);
[LM_sent_ar,cval_sent_ar]               = StabilityTestDyn(frec,[Lpmi(:,[1 2 4]) Lcc(:,[1 6])],trimr);
[LM_sent_dynar,cval_sent_dynar]         = StabilityTestDyn(frec,[Lpmi(:,[1 2 4]) Lcc(:,[1 6]) Llrec],trimr);

% Sentiment and classic predictor models
[LM_classent_stat,cval_classent_stat]   = StabilityTestStat(frec,[Lpmi(:,1) Lcc(:,1) Ltms(:,6) Lret(:,[2 4])],trimr);
[LM_classent_dyn,cval_classent_dyn]     = StabilityTestStat(frec,[Lpmi(:,1) Lcc(:,1) Ltms(:,6) Lret(:,[2 4]) Llrec],trimr);
[LM_classent_ar,cval_classent_ar]       = StabilityTestDyn(frec,[Lpmi(:,6) Lcc(:,[1 5]) Ltms(:,6) Lret(:,2)],trimr);
[LM_classent_dynar,cval_classent_dynar] = StabilityTestDyn(frec,[Lpmi(:,5) Lcc(:,[1 5]) Ltms(:,6) Lret(:,2) Llrec],trimr);

% Sentiment and factor predictor models
[LM_facsent_stat,cval_facsent_stat]     = StabilityTestStat(frec,[Lpmi(:,1) Lcc(:,1) Lfhat1(:,2) Lfhat6(:,6) Lfhat14(:,1)],trimr);
[LM_facsent_dyn,cval_facsent_dyn]       = StabilityTestStat(frec,[Lpmi(:,1) Lcc(:,1) Lfhat1(:,2) Lfhat6(:,6) Lfhat14(:,1) Llrec],trimr);
[LM_facsent_ar,cval_facsent_ar]         = StabilityTestDyn(frec,[Lpmi(:,1) Lfhat1(:,2) Lfhat2(:,2) Lfhat6(:,6) Lfhat14(:,1)],trimr);
[LM_facsent_dynar,cval_facsent_dynar]   = StabilityTestDyn(frec,[Lpmi(:,3) Lfhat1(:,2) Lfhat2(:,2) Lfhat6(:,6) Lfhat14(:,1) Llrec],trimr);


% Creating LaTeX Table for instability tests
fid = fopen('TeX/Tables/cem_stability_tests.tex','w');
fprintf(fid,'%s & %s & %s & %s \\\\\\midrule\n','','$\sup$-LM$_{t}\left(\omega\right)$','CV','H$_{0}$; No break');
fprintf(fid,'%s & %.2f & %.2f & %s \\\\\n','Model 1',LM_tms_stat,cval_tms_stat,'Not rejected');
fprintf(fid,'%s & %.2f & %.2f & %s \\\\\n','Model 2',LM_class_stat,cval_class_stat,'Not rejected');
fprintf(fid,'%s & %.2f & %.2f & %s \\\\\n','Model 3',LM_pmi_stat,cval_pmi_stat,'Not rejected');
fprintf(fid,'%s & %.2f & %.2f & %s \\\\\n','Model 4',LM_cc_stat,cval_cc_stat,'Not rejected');
fprintf(fid,'%s & %.2f & %.2f & %s \\\\\n','Model 5',LM_sent_stat,cval_sent_stat,'Not rejected');
fprintf(fid,'%s & %.2f & %.2f & %s \\\\\n','Model 6',LM_classent_stat,cval_classent_stat,'Not rejected');
fprintf(fid,'%s & %.2f & %.2f & %s \\\\\n','Model 7',LM_facsent_stat,cval_facsent_stat,'Not rejected');
fprintf(fid, '\n');
fclose(fid);

%-----------------------------------------------------------------------------------------------------------------------
%% CONDUCTING THE OUT-OF-SAMPLE FORECASTING EXERCISE
%-----------------------------------------------------------------------------------------------------------------------

disp('Conducting the out-of-sample forecasting exercise');

% Setting preliminaries
nhorz       = 5;
horz        = [1 3 6 9 12];
T           = size(frec,1);
R           = 235;                                      % In-Sample period 1978.6 - 1998.1 for explanatory variables
P           = T-R+1;                                    % Out-of-sample evaluation period

% Preallocations
nber_t                  = NaN(P,nhorz);
recprop_t               = NaN(P,1);
frcst_tms_stat          = NaN(P,nhorz);
frcst_pmi_stat          = NaN(P,nhorz);
frcst_cc_stat           = NaN(P,nhorz);
frcst_class_stat        = NaN(P,nhorz);
frcst_sent_stat         = NaN(P,nhorz);
frcst_classent_stat     = NaN(P,nhorz);
frcst_facsent_stat      = NaN(P,nhorz);
lps_tms_stat            = NaN(P,nhorz);
lps_pmi_stat            = NaN(P,nhorz);
lps_cc_stat             = NaN(P,nhorz);
lps_class_stat          = NaN(P,nhorz);
lps_sent_stat           = NaN(P,nhorz);
lps_classent_stat       = NaN(P,nhorz);
lps_facsent_stat        = NaN(P,nhorz);
lps_constant            = NaN(P,nhorz);

% Setting up script-based progress counter
reverseStr = '';
fprintf(1,'Progress of forecast procedure:  ');

% Out-of-sample forecasting loop
for hh = 1:nhorz
    
    h = horz(hh);
    
    for t=1:P-h
             
        msg = sprintf('Horizon %d Forecast %d', h, t);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        % Actual recession dates
        nber_t(t,hh) = frec(R+h+t-3,1);
        recprop_t(t,1) = mean(rec(1:R+t-2,1));

        % Term-spread predictor model
        rest = statProbit(frec(h:R+h+t-4,1),Ltms(1:R+t-3,6));
        frcst_tms_stat(t,hh) = norm_cdf([1 Ltms(R+t-1,6)]*rest.beta);
        
        % PMI predictor model
        rest = statProbit(frec(h:R+h+t-4,1),Lpmi(1:R+t-3,1));
        frcst_pmi_stat(t,hh) = norm_cdf([1 Lpmi(R+t-1,1)]*rest.beta);
  
        % CC predictor model
        rest = statProbit(frec(h:R+h+t-4,1),Lcc(1:R+t-3,1));
        frcst_cc_stat(t,hh) = norm_cdf([1 Lcc(R+t-1,1)]*rest.beta);
        
        % Classic predictor models
        rest = statProbit(frec(h:R+h+t-4,1),[Ltms(1:R+t-3,6) Lfed(1:R+t-3,6) Lret(1:R+t-3,[2 4 6])]);
        frcst_class_stat(t,hh) = norm_cdf([1 Ltms(R+t-1,6) Lfed(R+t-1,6) Lret(R+t-1,[2 4 6])]*rest.beta);

        % Sentiment predictor models
        rest = statProbit(frec(h:R+h+t-4,1),[Lpmi(1:R+t-3,1) Lcc(1:R+t-3,[1 6])]);
        frcst_sent_stat(t,hh) = norm_cdf([1 Lpmi(R+t-1,1) Lcc(R+t-1,[1 6])]*rest.beta);
        
        % Sentiment and classic predictor models
        rest = statProbit(frec(h:R+h+t-4,1),[Lpmi(1:R+t-3,1) Lcc(1:R+t-3,1) Ltms(1:R+t-3,6) Lret(1:R+t-3,[2 4])]);
        frcst_classent_stat(t,hh) = norm_cdf([1 Lpmi(R+t-1,1) Lcc(R+t-1,1) Ltms(R+t-1,6) Lret(R+t-1,[2 4])]*rest.beta);
        
        % Sentiment and factor predictor models
        [Lfhat1t,Lfhat2t,Lfhat3t,Lfhat4t,Lfhat5t,Lfhat6t,Lfhat7t,Lfhat8t,Lfhat9t,Lfhat10t,Lfhat11t,Lfhat12t,Lfhat13t,Lfhat14t,Lfhat15t]...
        = recursive_factor_construction(rawdata(1:R+t+7,:),tcode,15);
    
        rest = statProbit(frec(h:R+h+t-4,1),[Lpmi(1:R+t-3,1) Lcc(1:R+t-3,1) Lfhat1(1:R+t-3,2) Lfhat6(1:R+t-3,6) Lfhat14(1:R+t-3,1)]);
        frcst_facsent_stat(t,hh) = norm_cdf([1 Lpmi(R+t-1,1) Lcc(R+t-1,1) Lfhat1(R+t-1,2) Lfhat6(R+t-1,6) Lfhat14(R+t-1,1)]*rest.beta);

    end
      
    % Correcting for zeros and ones
    tolz = 1e-6; 
    frcst_tms_stat(:,hh)        = min(max(frcst_tms_stat(:,hh),tolz),1-tolz);
    frcst_pmi_stat(:,hh)        = min(max(frcst_pmi_stat(:,hh),tolz),1-tolz);
    frcst_cc_stat(:,hh)         = min(max(frcst_cc_stat(:,hh),tolz),1-tolz);
    frcst_class_stat(:,hh)      = min(max(frcst_class_stat(:,hh),tolz),1-tolz);
    frcst_sent_stat(:,hh)       = min(max(frcst_sent_stat(:,hh),tolz),1-tolz);
    frcst_classent_stat(:,hh)   = min(max(frcst_classent_stat(:,hh),tolz),1-tolz);
    frcst_facsent_stat(:,hh)    = min(max(frcst_facsent_stat(:,hh),tolz),1-tolz);

    % Computing QPS (MSE)
    qps_tms_stat(1,hh)          = 2*mean((frcst_tms_stat(1:end-h,hh) - nber_t(1:end-h,hh)).^2);
    qps_pmi_stat(1,hh)          = 2*mean((frcst_pmi_stat(1:end-h,hh) - nber_t(1:end-h,hh)).^2);
    qps_cc_stat(1,hh)           = 2*mean((frcst_cc_stat(1:end-h,hh) - nber_t(1:end-h,hh)).^2);
    qps_class_stat(1,hh)        = 2*mean((frcst_class_stat(1:end-h,hh) - nber_t(1:end-h,hh)).^2);
    qps_sent_stat(1,hh)         = 2*mean((frcst_sent_stat(1:end-h,hh) - nber_t(1:end-h,hh)).^2);
    qps_classent_stat(1,hh)     = 2*mean((frcst_classent_stat(1:end-h,hh) - nber_t(1:end-h,hh)).^2);
    qps_facsent_stat(1,hh)      = 2*mean((frcst_facsent_stat(1:end-h,hh) - nber_t(1:end-h,hh)).^2);

    % Computing LPS
    lps_constant(1:P-h,hh)      =  -(nber_t(1:end-h,hh).*log(recprop_t(1:end-h,1)) + (1-nber_t(1:end-h,hh)).*log(1-recprop_t(1:end-h,1)));
    lps_tms_stat(1:P-h,hh)      = -(nber_t(1:end-h,hh).*log(frcst_tms_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_tms_stat(1:end-h,hh)));
    lps_pmi_stat(1:P-h,hh)      = -(nber_t(1:end-h,hh).*log(frcst_pmi_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_pmi_stat(1:end-h,hh)));
    lps_cc_stat(1:P-h,hh)       = -(nber_t(1:end-h,hh).*log(frcst_cc_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_cc_stat(1:end-h,hh)));
    lps_class_stat(1:P-h,hh)    = -(nber_t(1:end-h,hh).*log(frcst_class_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_class_stat(1:end-h,hh)));
    lps_sent_stat(1:P-h,hh)     = -(nber_t(1:end-h,hh).*log(frcst_sent_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_sent_stat(1:end-h,hh)));
    lps_classent_stat(1:P-h,hh) = -(nber_t(1:end-h,hh).*log(frcst_classent_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_classent_stat(1:end-h,hh)));
    lps_facsent_stat(1:P-h,hh)  = -(nber_t(1:end-h,hh).*log(frcst_facsent_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_facsent_stat(1:end-h,hh)));

    % Computing pseudo-R2
    LR  = sum(nber_t(1:end-h,hh).*log(recprop_t(1:end-h,1)) + (1-nber_t(1:end-h,hh)).*log(1-recprop_t(1:end-h,1)));
    LU  = sum(nber_t(1:end-h,hh).*log(frcst_tms_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_tms_stat(1:end-h,hh)));
    r2_tms_stat(1,hh) = 1-(LU/LR)^((-2/(P-h))*LR);
    LU = sum(nber_t(1:end-h,hh).*log(frcst_pmi_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_pmi_stat(1:end-h,hh)));
    r2_pmi_stat(1,hh) = 1-(LU/LR)^((-2/(P-h))*LR);
    LU = sum(nber_t(1:end-h,hh).*log(frcst_cc_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_cc_stat(1:end-h,hh)));
    r2_cc_stat(1,hh) = 1-(LU/LR)^((-2/(P-h))*LR);
    LU = sum(nber_t(1:end-h,hh).*log(frcst_class_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_class_stat(1:end-h,hh)));
    r2_class_stat(1,hh) = 1-(LU/LR)^((-2/(P-h))*LR);
    LU = sum(nber_t(1:end-h,hh).*log(frcst_sent_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_sent_stat(1:end-h,hh)));
    r2_sent_stat(1,hh) = 1-(LU/LR)^((-2/(P-h))*LR);
    LU = sum(nber_t(1:end-h,hh).*log(frcst_classent_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_classent_stat(1:end-h,hh)));
    r2_classent_stat(1,hh) = 1-(LU/LR)^((-2/(P-h))*LR);
    LU = sum(nber_t(1:end-h,hh).*log(frcst_facsent_stat(1:end-h,hh)) + (1-nber_t(1:end-h,hh)).*log(1-frcst_facsent_stat(1:end-h,hh)));
    r2_facsent_stat(1,hh) = 1-(LU/LR)^((-2/(P-h))*LR);
    clear LU

    % Computing correct recession classification
    indx = find(nber_t(1:end-h,hh)==1);
    cr50_tms_stat(1,hh)         = sum(((frcst_tms_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr50_pmi_stat(1,hh)         = sum(((frcst_pmi_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr50_cc_stat(1,hh)          = sum(((frcst_cc_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr50_class_stat(1,hh)       = sum(((frcst_class_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr50_sent_stat(1,hh)        = sum(((frcst_sent_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr50_classent_stat(1,hh)    = sum(((frcst_classent_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr50_facsent_stat(1,hh)     = sum(((frcst_facsent_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 

    cr25_tms_stat(1,hh)         = sum(((frcst_tms_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr25_pmi_stat(1,hh)         = sum(((frcst_pmi_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr25_cc_stat(1,hh)          = sum(((frcst_cc_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr25_class_stat(1,hh)       = sum(((frcst_class_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr25_sent_stat(1,hh)        = sum(((frcst_sent_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr25_classent_stat(1,hh)    = sum(((frcst_classent_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    cr25_facsent_stat(1,hh)     = sum(((frcst_facsent_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 

    % Computing correct expansion classification
    indx = find(nber_t(1:end-h,hh)==0);
    ce50_tms_stat(1,hh)         = sum(((frcst_tms_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce50_pmi_stat(1,hh)         = sum(((frcst_pmi_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce50_cc_stat(1,hh)          = sum(((frcst_cc_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce50_class_stat(1,hh)       = sum(((frcst_class_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce50_sent_stat(1,hh)        = sum(((frcst_sent_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce50_classent_stat(1,hh)    = sum(((frcst_classent_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1);  
    ce50_facsent_stat(1,hh)     = sum(((frcst_facsent_stat(indx,hh) > 0.5) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 

    ce25_tms_stat(1,hh)         = sum(((frcst_tms_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce25_pmi_stat(1,hh)         = sum(((frcst_pmi_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce25_cc_stat(1,hh)          = sum(((frcst_cc_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce25_class_stat(1,hh)       = sum(((frcst_class_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce25_sent_stat(1,hh)        = sum(((frcst_sent_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce25_classent_stat(1,hh)    = sum(((frcst_classent_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 
    ce25_facsent_stat(1,hh)     = sum(((frcst_facsent_stat(indx,hh) > 0.25) == nber_t(indx,hh)))/size(nber_t(indx,hh),1); 

    % Computing AUROC
    [~,~,~,auc_tms_stat(1,hh)]      = perfcurve(nber_t(1:end-h,hh),frcst_tms_stat(1:end-h,hh),1);
    [~,~,~,auc_pmi_stat(1,hh)]      = perfcurve(nber_t(1:end-h,hh),frcst_pmi_stat(1:end-h,hh),1);
    [~,~,~,auc_cc_stat(1,hh)]       = perfcurve(nber_t(1:end-h,hh),frcst_cc_stat(1:end-h,hh),1);
    [~,~,~,auc_class_stat(1,hh)]    = perfcurve(nber_t(1:end-h,hh),frcst_class_stat(1:end-h,hh),1);
    [~,~,~,auc_sent_stat(1,hh)]     = perfcurve(nber_t(1:end-h,hh),frcst_sent_stat(1:end-h,hh),1);
    [~,~,~,auc_classent_stat(1,hh)] = perfcurve(nber_t(1:end-h,hh),frcst_classent_stat(1:end-h,hh),1);
    [~,~,~,auc_facsent_stat(1,hh)]  = perfcurve(nber_t(1:end-h,hh),frcst_facsent_stat(1:end-h,hh),1);

    % Computing standard errors and testing AUC null hypothesis
    n0          = size(nber_t(indx,hh),1);
    n1          = P-h-n0;
    Q1          = auc_tms_stat(1,hh) / (2-auc_tms_stat(1,hh));
    Q2          = (2*auc_tms_stat(1,hh)^2) / (1+auc_tms_stat(1,hh)); 
    AUCsigma    = sqrt((1/(n0*n1))*(auc_tms_stat(1,hh)*(1-auc_tms_stat(1,hh))+(n1-1)*...
    (Q1-auc_tms_stat(1,hh)^2)+(n0-1)*(Q2-auc_tms_stat(1,hh)^2)));
    aucp_tms_stat(1,hh) = (1-normcdf((auc_tms_stat(1,hh)-0.5)/AUCsigma,0,1));

    Q1          = auc_pmi_stat(1,hh) / (2-auc_pmi_stat(1,hh));
    Q2          = (2*auc_pmi_stat(1,hh)^2) / (1+auc_pmi_stat(1,hh)); 
    AUCsigma    = sqrt((1/(n0*n1))*(auc_pmi_stat(1,hh)*(1-auc_pmi_stat(1,hh))+(n1-1)*...
    (Q1-auc_pmi_stat(1,hh)^2)+(n0-1)*(Q2-auc_pmi_stat(1,hh)^2)));
    aucp_pmi_stat(1,hh) = (1-normcdf((auc_pmi_stat(1,hh)-0.5)/AUCsigma,0,1));

    Q1          = auc_cc_stat(1,hh) / (2-auc_cc_stat(1,hh));
    Q2          = (2*auc_cc_stat(1,hh)^2) / (1+auc_cc_stat(1,hh)); 
    AUCsigma    = sqrt((1/(n0*n1))*(auc_cc_stat(1,hh)*(1-auc_cc_stat(1,hh))+(n1-1)*...
    (Q1-auc_cc_stat(1,hh)^2)+(n0-1)*(Q2-auc_cc_stat(1,hh)^2)));
    aucp_cc_stat(1,hh) = (1-normcdf((auc_cc_stat(1,hh)-0.5)/AUCsigma,0,1));

    Q1          = auc_class_stat(1,hh) / (2-auc_class_stat(1,hh));
    Q2          = (2*auc_class_stat(1,hh)^2) / (1+auc_class_stat(1,hh)); 
    AUCsigma    = sqrt((1/(n0*n1))*(auc_class_stat(1,hh)*(1-auc_class_stat(1,hh))+(n1-1)*...
    (Q1-auc_class_stat(1,hh)^2)+(n0-1)*(Q2-auc_class_stat(1,hh)^2)));
    aucp_class_stat(1,hh) = (1-normcdf((auc_class_stat(1,hh)-0.5)/AUCsigma,0,1));

    Q1          = auc_sent_stat(1,hh) / (2-auc_sent_stat(1,hh));
    Q2          = (2*auc_sent_stat(1,hh)^2) / (1+auc_sent_stat(1,hh)); 
    AUCsigma    = sqrt((1/(n0*n1))*(auc_sent_stat(1,hh)*(1-auc_sent_stat(1,hh))+(n1-1)*...
    (Q1-auc_sent_stat(1,hh)^2)+(n0-1)*(Q2-auc_sent_stat(1,hh)^2)));
    aucp_sent_stat(1,hh) = (1-normcdf((auc_sent_stat(1,hh)-0.5)/AUCsigma,0,1));

    Q1          = auc_classent_stat(1,hh) / (2-auc_classent_stat(1,hh));
    Q2          = (2*auc_classent_stat(1,hh)^2) / (1+auc_classent_stat(1,hh)); 
    AUCsigma    = sqrt((1/(n0*n1))*(auc_classent_stat(1,hh)*(1-auc_classent_stat(1,hh))+(n1-1)*...
    (Q1-auc_classent_stat(1,hh)^2)+(n0-1)*(Q2-auc_classent_stat(1,hh)^2)));
    aucp_classent_stat(1,hh) = (1-normcdf((auc_classent_stat(1,hh)-0.5)/AUCsigma,0,1));

    Q1          = auc_facsent_stat(1,hh) / (2-auc_facsent_stat(1,hh));
    Q2          = (2*auc_facsent_stat(1,hh)^2) / (1+auc_facsent_stat(1,hh)); 
    AUCsigma    = sqrt((1/(n0*n1))*(auc_facsent_stat(1,hh)*(1-auc_facsent_stat(1,hh))+(n1-1)*...
    (Q1-auc_facsent_stat(1,hh)^2)+(n0-1)*(Q2-auc_facsent_stat(1,hh)^2)));
    aucp_facsent_stat(1,hh) = (1-normcdf((auc_facsent_stat(1,hh)-0.5)/AUCsigma,0,1));


    % Computing Diebold-Mariano tests
    [DM_class_classent(1,hh),DMp_class_classent(1,hh)]      = perform_dm_lps_test(lps_class_stat(1:P-h,hh),lps_classent_stat(1:P-h,hh));
    [DM_class_sent(1,hh),DMp_class_sent(1,hh)]              = perform_dm_lps_test(lps_class_stat(1:P-h,hh),lps_sent_stat(1:P-h,hh));
    [DM_class_facsent(1,hh),DMp_class_facsent(1,hh)]        = perform_dm_lps_test(lps_class_stat(1:P-h,hh),lps_facsent_stat(1:P-h,hh));
    [DM_sent_facsent(1,hh),DMp_sent_facsent(1,hh)]          = perform_dm_lps_test(lps_sent_stat(1:P-h,hh),lps_facsent_stat(1:P-h,hh));
    [DM_classsent_facsent(1,hh),DMp_classent_facsent(1,hh)] = perform_dm_lps_test(lps_classent_stat(1:P-h,hh),lps_facsent_stat(1:P-h,hh));

    % Computing test for difference in AUC
    [dauc_class_classent(1,hh),daucp_class_classent(1,hh)]      = perform_auc_test(nber_t(1:end-h,hh),auc_classent_stat(hh),auc_class_stat(hh),frcst_classent_stat(1:end-h,hh),frcst_class_stat(1:end-h,hh));
    [dauc_class_sent(1,hh),daucp_class_sent(1,hh)]              = perform_auc_test(nber_t(1:end-h,hh),auc_sent_stat(hh),auc_class_stat(hh),frcst_sent_stat(1:end-h,hh),frcst_class_stat(1:end-h,hh));
    [dauc_class_facsent(1,hh),daucp_class_facsent(1,hh)]        = perform_auc_test(nber_t(1:end-h,hh),auc_facsent_stat(hh),auc_class_stat(hh),frcst_facsent_stat(1:end-h,hh),frcst_class_stat(1:end-h,hh));
    [dauc_sent_facsent(1,hh),daucp_sent_facsent(1,hh)]          = perform_auc_test(nber_t(1:end-h,hh),auc_facsent_stat(hh),auc_sent_stat(hh),frcst_facsent_stat(1:end-h,hh),frcst_sent_stat(1:end-h,hh));
    [dauc_classent_facsent(1,hh),daucp_classent_facsent(1,hh)]  = perform_auc_test(nber_t(1:end-h,hh),auc_facsent_stat(hh),auc_classent_stat(hh),frcst_facsent_stat(1:end-h,hh),frcst_classent_stat(1:end-h,hh));
      
end
fprintf('\n');

% Creating LaTeX Table for out-of-sample results
fid = fopen('TeX/Tables/cem_out_of_sample_results.tex','w');
fprintf(fid,'%s & %.0f & %.0f & %.0f & %.0f & %.0f & %s & %.0f & %.0f & %.0f & %.0f & %.0f \\\\\\midrule\n','',[1 3 6 9 12],'',[1 3 6 9 12]);
fprintf(fid,'%s & %s & %s & %s \\\\\\cmidrule{2-6}\\cmidrule{8-12}\n','','\multicolumn{5}{c}{\textit{pseudo-R$^{2}$}}','','\multicolumn{5}{c}{\textit{LPS}}');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 1',r2_tms_stat,'',nanmean(lps_tms_stat));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 2',r2_class_stat,'',nanmean(lps_class_stat));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 3',r2_pmi_stat,'',nanmean(lps_pmi_stat));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 4',r2_cc_stat,'',nanmean(lps_cc_stat));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 5',r2_sent_stat,'',nanmean(lps_sent_stat));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 6',r2_classent_stat,'',nanmean(lps_classent_stat));
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\\cmidrule{2-6}\\cmidrule{8-12}\n\n','Model 7',r2_facsent_stat,'',nanmean(lps_facsent_stat));
fprintf(fid,'%s & %s & %s & %s \\\\\\cmidrule{2-6}\\cmidrule{8-12}\n','','\multicolumn{5}{c}{\textit{AUC}}','','\multicolumn{5}{c}{\textit{p-val}}');
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 1',auc_tms_stat,'',aucp_tms_stat);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 2',auc_class_stat,'',aucp_class_stat);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 3',auc_pmi_stat,'',aucp_pmi_stat);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 4',auc_cc_stat,'',aucp_cc_stat);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 5',auc_sent_stat,'',aucp_sent_stat);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 6',auc_classent_stat,'',aucp_classent_stat);
fprintf(fid,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n','Model 7',auc_facsent_stat,'',aucp_facsent_stat);
fprintf(fid,'\n');
fclose(fid);

%-----------------------------------------------------------------------------------------------------------------------
%% PLOTTING OUT-OF-SAMPLE FIT OF SELECTED PROBIT MODELS
%-----------------------------------------------------------------------------------------------------------------------

disp('Plotting out-of-sample fit of selected probit models');

% Lopping over forecast horizons
for hh = 1:5

    % Set forecast horizon
    h = horz(hh);

    % Collect model forecasts
    prob_rec = [
        frcst_class_stat(1:end-h,hh) ...
        frcst_sent_stat(1:end-h,hh) ...
        frcst_classent_stat(1:end-h,hh)...
        frcst_facsent_stat(1:end-h,hh)
    ];

    % Plot model implied recession probabilities
    figure;
    for iFig = 1:4

        subplot(2,2,iFig);
        hold on

        b1 = bar(datenum_indx(241:end-h,:),nber_t(1:end-h,hh));
        b1.EdgeColor = [0.8 0.8 0.8];
        b1.FaceColor = [0.8 0.8 0.8];
        b1.BarWidth  = 1;

        p1 = plot(datenum_indx(241:end-h,:),prob_rec(:,iFig));
        p1.LineStyle = '-'; p1.Color = 'k';
        thinl = line([datenum(1998,1,1) ; datenum(2012,1,1)],[0.5 ; 0.5],'LineStyle',':','Color','k');
        set(get(get(thinl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        datetick('x','yyyy');
        axis([-inf inf 0 1.004]);
        title(model_str{iFig},'FontSize',8);
        xlabel('Time');
        ylabel('Probability');
        box on

    end
    set(gcf,'PaperUnits','Centimeters','PaperSize',[21 15],'PaperPosition',[0.5 0.5 20 14]);
    print(gcf,'-depsc',['TeX/Figures/Figure3_' num2str(h) '.eps']);

end

%-----------------------------------------------------------------------------------------------------------------------
%% PLOTTING RECEIVER OPERATING CHARACTERISTICS (ROC) CURVES
%-----------------------------------------------------------------------------------------------------------------------

disp('Plotting receiver operating characteristics (ROC) curves');

% Plot in-sample ROC curves
figure;
hold on

[X,Y] = perfcurve(frec,res_class_stat.yhat,1);
p1 = plot(X,Y);
p1.LineStyle = ':';
p1.Color = 'k';
p1.LineWidth = 0.8;
p1.MarkerSize = 4;
hs = line('XData', [0 1], 'YData', [0 1], 'LineStyle', '-','LineWidth', 1, 'Color','k');
set(get(get(hs, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

[X,Y] = perfcurve(frec,res_sent_stat.yhat,1);
p1 = plot(X,Y);
p1.LineStyle = '-.';
p1.Color = 'k';
p1.LineWidth = 0.8;
p1.MarkerSize = 4;

[X,Y] = perfcurve(frec,res_classent_stat.yhat,1);
p1 = plot(X,Y);
p1.LineStyle = '--';
p1.Color = 'k';
p1.LineWidth = 0.8;
p1.MarkerSize = 4;

[X,Y] = perfcurve(frec,res_facsent_stat.yhat,1);
p1 = plot(X,Y);
p1.LineStyle = '-';
p1.Color = 'k';
p1.LineWidth = 0.8;
p1.MarkerSize = 4;
hold off
box on

xlabel('False positive rate'); 
ylabel('True positive rate')
leg = legend(model_str,'Location','SouthEast');
set(leg,'box','off');
print(gcf,'-depsc','TeX/Figures/Figure4.eps');

% Plotting out-of-sample ROC curves
for hh = 1:nhorz
    
    % Set forecast horizon
    h = horz(hh);

    % Creating figures
    figure;
    hold on

    [X,Y] = perfcurve(nber_t(1:end-h,hh),frcst_class_stat(1:end-h,hh),1);
    p1 = plot(X,Y);
    p1.LineStyle = ':';
    p1.Color = 'k';
    p1.LineWidth = 0.8;
    p1.MarkerSize = 4;
    hs = line('XData', [0 1], 'YData', [0 1], 'LineStyle', '-','LineWidth', 1, 'Color','k');
    set(get(get(hs, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    [X,Y] = perfcurve(nber_t(1:end-h,hh),frcst_sent_stat(1:end-h,hh),1);
    p1 = plot(X,Y);
    p1.LineStyle = '-.';
    p1.Color = 'k';
    p1.LineWidth = 0.8;
    p1.MarkerSize = 4;

    [X,Y] = perfcurve(nber_t(1:end-h,hh),frcst_classent_stat(1:end-h,hh),1);
    p1 = plot(X,Y);
    p1.LineStyle = '--';
    p1.Color = 'k';
    p1.LineWidth = 0.8;
    p1.MarkerSize = 4;

    [X,Y] = perfcurve(nber_t(1:end-h,hh),frcst_facsent_stat(1:end-h,hh),1);
    p1 = plot(X,Y);
    p1.LineStyle = '-';
    p1.Color = 'k';
    p1.LineWidth = 0.8;
    p1.MarkerSize = 4;
    hold off
    box on

    xlabel('False positive rate'); 
    ylabel('True positive rate');
    legend('Model 2','Model 5','Model 6','Model 7','Location','SouthEast');
    print(gcf,'-depsc',['TeX/Figures/Figure5_' num2str(h) '.eps']);

end

%-----------------------------------------------------------------------------------------------------------------------
%% ESTIMATING LINEAR REGRESSION MODELS
%-----------------------------------------------------------------------------------------------------------------------

disp('Estimating linear regression models for industrial production growth');

% Setting up the data
mdata = [
    ip_growth ...
    Lpmi(:,1) ...
    Lpmi(:,3) ...
    Lcc(:,1) ...
    Ltms(:,4) ...
    Lfed(:,4) ...
    Lret(:,2) ...
    Lret(:,3) ...
    Lret(:,6) ...
    Lfhat1(:,3) ...
    Lfhat6(:,6)
];

% Estimating predictive regression models
res_tms_stat        = nwest(ip_growth,[ones(402,1) Ltms(:,4)],3);
res_pmi_stat        = nwest(ip_growth,[ones(402,1) Lpmi(:,1)],3);
res_cc_stat         = nwest(ip_growth,[ones(402,1) Lcc(:,1)],3);
res_class_stat      = nwest(ip_growth,[ones(402,1) Ltms(:,4) Lfed(:,4) Lret(:,[2 3 6])],3);
res_sent_stat       = nwest(ip_growth,[ones(402,1) Lpmi(:,[1 3]) Lcc(:,1)],3);
res_classent_stat   = nwest(ip_growth,[ones(402,1) Lpmi(:,[1 3]) Lcc(:,1) Lret(:,[2 6])],3);
res_facsent_stat    = nwest(ip_growth,[ones(402,1) Lpmi(:,[1 3]) Lcc(:,1) Lfhat1(:,3) Lfhat6(:,6)],3);

% Create latex table
strvar = {'PMI$_{t-1}$','PMI$_{t-3}$','CC$_{t-1}$','TS$_{t-4}$','FFR$_{t-4}$','RET$_{t-2}$','RET$_{t-3}$','RET$_{t-6}$','$\hat{f}_{1,t-3}$','$\hat{f}_{6,t-6}$'};
FID = fopen('TeX/Tables/cem_forecasting_industrial_production.tex','w');
fprintf(FID,'%s & %s & %s & %s & %s & %s & %s & %s \\\\ \\midrule \n','Variable','Model 1','Model 2','Model 3','Model 4','Model 5','Model 6','Model 7');
fprintf(FID,'%s & %s & %s & %.2f & %s & %.2f & %.2f & %.2f \\\\',strvar{1},'','',res_pmi_stat.beta(2),'',res_sent_stat.beta(2),res_classent_stat.beta(2),res_facsent_stat.beta(2));
fprintf(FID,'%s & %s & %s & (%.2f) & %s & (%.2f) & (%.2f) & (%.2f) \\\\','','','',res_pmi_stat.stdb(2),'',res_sent_stat.stdb(2),res_classent_stat.stdb(2),res_facsent_stat.stdb(2));
fprintf(FID,'%s & %s & %s & %s & %s & %.2f & %.2f & %.2f \\\\',strvar{2},'','','','',res_sent_stat.beta(3),res_classent_stat.beta(3),res_facsent_stat.beta(3));
fprintf(FID,'%s & %s & %s & %s & %s & (%.2f) & (%.2f) & (%.2f) \\\\','','','','','',res_sent_stat.stdb(3),res_classent_stat.stdb(3),res_facsent_stat.stdb(3));
fprintf(FID,'%s & %s & %s & %s & %.2f & %.2f & %.2f & %.2f \\\\',strvar{3},'','','',res_cc_stat.beta(2),res_sent_stat.beta(4),res_classent_stat.beta(4),res_facsent_stat.beta(4));
fprintf(FID,'%s & %s & %s & %s & (%.2f) & (%.2f) & (%.2f) & (%.2f) \\\\','','','','',res_cc_stat.stdb(2),res_sent_stat.stdb(4),res_classent_stat.stdb(4),res_facsent_stat.stdb(4));
fprintf(FID,'%s & %.2f & %.2f & %s & %s & %s & %s & %s \\\\',strvar{4},res_tms_stat.beta(2),res_class_stat.beta(2),'','','','','');
fprintf(FID,'%s & (%.2f) & (%.2f) & %s & %s & %s & %s & %s \\\\','',res_tms_stat.stdb(2),res_class_stat.stdb(2),'','','','','');
fprintf(FID,'%s & %s & %.2f & %s & %s & %s & %s & %s \\\\',strvar{5},'',res_class_stat.beta(3),'','','','','');
fprintf(FID,'%s & %s & (%.2f) & %s & %s & %s & %s & %s \\\\','','',res_class_stat.stdb(3),'','','','','');
fprintf(FID,'%s & %s & %.2f & %s & %s & %s & %.2f & %s \\\\',strvar{6},'',res_class_stat.beta(4),'','','',res_classent_stat.beta(5),'');
fprintf(FID,'%s & %s & (%.2f) & %s & %s & %s & (%.2f) & %s \\\\','','',res_class_stat.stdb(4),'','','',res_classent_stat.stdb(5),'');
fprintf(FID,'%s & %s & %.2f & %s & %s & %s & %s & %s \\\\',strvar{7},'',res_class_stat.beta(5),'','','','','');
fprintf(FID,'%s & %s & (%.2f) & %s & %s & %s & %s & %s \\\\','','',res_class_stat.stdb(5),'','','','','');
fprintf(FID,'%s & %s & %.2f & %s & %s & %s & %.2f & %s \\\\',strvar{8},'',res_class_stat.beta(6),'','','',res_classent_stat.beta(6),'');
fprintf(FID,'%s & %s & (%.2f) & %s & %s & %s & (%.2f) & %s \\\\','','',res_class_stat.stdb(6),'','','',res_classent_stat.stdb(6),'');
fprintf(FID,'%s & %s & %s & %s & %s & %s & %s & %.2f \\\\',strvar{9},'','','','','','',res_facsent_stat.beta(5));
fprintf(FID,'%s & %s & %s & %s & %s & %s & %s & (%.2f) \\\\','','','','','','','',res_facsent_stat.stdb(5));
fprintf(FID,'%s & %s & %s & %s & %s & %s & %s & %.2f \\\\',strvar{10},'','','','','','',res_facsent_stat.beta(6));
fprintf(FID,'%s & %s & %s & %s & %s & %s & %s & (%.2f) \\\\','','','','','','','',res_facsent_stat.stdb(6));
fprintf(FID,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\','BIC',res_tms_stat.sic,res_class_stat.sic,res_pmi_stat.sic,res_cc_stat.sic,res_sent_stat.sic,res_classent_stat.sic,res_facsent_stat.sic);
fprintf(FID,'%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\','adj. R$^{2}$',res_tms_stat.rbar,res_class_stat.rbar,res_pmi_stat.rbar,res_cc_stat.rbar,res_sent_stat.rbar,res_classent_stat.rbar,res_facsent_stat.rbar);
fprintf(FID, '\n');
fclose(FID);

%-----------------------------------------------------------------------------------------------------------------------
%% COMPUTING CODE RUN TIME
%-----------------------------------------------------------------------------------------------------------------------

tEnd = toc(tStart);
fprintf('Runtime: %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
disp('Routine Completed');

%-----------------------------------------------------------------------------------------------------------------------
% END OF SCRIPT
%-----------------------------------------------------------------------------------------------------------------------