## Forecasting US Recessions: The Role of Sentiment

* Author: Charlotte Christiansen, Jonas Nygaard Eriksen, and Stig Vinther Møller
* Maintainer: Jonas Nygaard Eriksen (jeriksen (at) econ.au.dk)
* Link to paper: [Forecasting US Recessions: The Role of Sentiment](http://dx.doi.org/10.1016/j.jbankfin.2014.06.017)
* Paper status: Published in Journal of Banking and Finance

### Abstract
> We study the role of sentiment variables as predictors for US recessions. We combine sentiment variables with either classical recession predictors or common factors based on a large panel of macroeconomic and financial variables. Sentiment variables hold vast predictive power for US recessions in excess of both the classical recession predictors and the common factors. The strong importance of the sentiment vari- ables is documented both in-sample and out-of-sample.

### Replication files
This repository contains the source files necessary for replicating the empirical results presented in the paper *_Forecasting US recessions: The role of sentiment_*. Specifically, this repository contains all Matlab files used for generating the empirical results as well as the tables and figures appearing in the paper. 

#### Generate figures and tables
To re-create the results, simply run `main_cem.m` in Matlab. It calls all dependencies and stores `TeX` tables and `eps` figures in the `TeX/Figures` and `TeX/Tables` folders. The data is stored in `.mat` files in the `Data/` folder, which contains further notes on the data. The code is run using Matlab R2015b for Mac. 

#### Compiling the LaTeX file
Open the `figures_tables.tex` file in the `TeX/` folder in our favorite editor and hit compile. The document is fairly plain vanilla, so hopefully there should be no trouble compiling. It loads figures and tables directly from the `TeX/Figures` and `TeX/Tables` folders, where the Matlab script `main_brp.m` stores the results. The LaTeX document is written and prepared in [Sublime Text 3](http://www.sublimetext.com/3) using the [LaTeX Tools](https://github.com/SublimeText/LaTeXTools) plug-in. 
