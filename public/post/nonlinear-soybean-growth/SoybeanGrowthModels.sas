*CLEARS SAS LOG AND RESULTS FOR CLEANER WORKING ENVIRONMENT;
dm "log; clear; odsresults; clear;";

* ----------------------------------------------------------------------------------------------------------------------;
* IMPORT ORIGINAL DATA ---------------------------------------------------------------------------------------;
* ----------------------------------------------------------------------------------------------------------------------;
PROC IMPORT
	DATAFILE = 'C:\Users\EmilyARobinson\Dropbox\Nonlinear\Soybean Growth\Data\soybean_data.csv'
	OUT = soybean_data
	REPLACE;
RUN;

TITLE "Soybean Data";
PROC PRINT DATA = soybean_data (OBS = 10) NOOBS;
RUN;

PROC SGPLOT DATA = soybean_data;
	SCATTER x = Days y = Leaf_Weight / group = Genotype;
RUN;

* ----------------------------------------------------------------------------------------------------------------------;
* 3-PARAMETER LOGISTIC OLS MODEL -------------------------------------------------------------------------;
* ----------------------------------------------------------------------------------------------------------------------;
TITLE "3- Parameter Logistic OLS Model";
*ODS SELECT ParameterEstimates;
PROC NLIN DATA = soybean_data;
	PARMS a = 20, b = 700, c = 0.125;
	MODEL Leaf_Weight = a/(1+b*exp(-c*Days));
	OUTPUT OUT = ypred P = yhat;
RUN;

PROC SORT DATA = ypred;
	BY Days;
RUN;

PROC SGPLOT DATA = ypred;
	SCATTER x = Days y = Leaf_Weight / group = Genotype;
	SERIES    x = Days y = yhat;
RUN;

* ----------------------------------------------------------------------------------------------------------------------;
* 4-PARAMETER LOGISTIC OLS MODEL -------------------------------------------------------------------------;
* ----------------------------------------------------------------------------------------------------------------------;
TITLE "4- Parameter Logistic OLS Model";
*ODS SELECT ParameterEstimates;
PROC NLIN DATA = soybean_data METHOD = marquardt;
	PARMS a = 0.2 b = 20 c = 50 d = 8;
	MODEL Leaf_Weight = a + (b-a)/(1+exp((c-Days)/d));
	OUTPUT OUT = ypred P = yhat;
RUN;

PROC SORT DATA = ypred;
	BY Days;
RUN;

PROC SGPLOT DATA = ypred;
	SCATTER x = Days y = Leaf_Weight / group = Genotype;
	SERIES    x = Days y = yhat;
RUN;

* ----------------------------------------------------------------------------------------------------------------------;
* 3-PARAMETER LOGISTIC OLS MODEL : BY GENOTYPE -----------------------------------------------------;
* ----------------------------------------------------------------------------------------------------------------------;
PROC IMPORT
	DATAFILE = 'C:\Users\EmilyARobinson\Dropbox\Nonlinear\Soybean Growth\Data\soybean_data2.csv'
	OUT = soybean_data2
	REPLACE;
RUN;

TITLE "Soybean Data - Indicator Column";
PROC PRINT DATA = soybean_data2 (OBS = 10) NOOBS;
RUN;

* BY GENOTYPE------------------------------------------------;
TITLE "3- Parameter Logistic OLS Model: By Genotype";
*ODS SELECT ParameterEstimates;
PROC NLIN DATA =soybean_data2 METHOD = marquardt;
	BY Genotype;
	PARMS a = 20 b = 700 c = 0.125;
	MODEL Leaf_Weight = a/(1+b*exp(-c*Days));
	OUTPUT OUT = ypred P = yhat;
RUN;
PROC SORT DATA = ypred;
	BY Genotype Days;
RUN;
PROC SGPLOT DATA = ypred;
	SCATTER x = Days y = Leaf_Weight / group = Genotype;
	SERIES    x = Days y = yhat / group = Genotype;
RUN;

* FULL MODEL ----------------------------------------------------;
PROC NLIN DATA = soybean_data2;
	PARMS a = 16, a1 = 4.78 b = 1035, b1 = -490 c = 0.125 c1 = -0.01;
	MODEL Leaf_Weight = (a + a1*Genotype_P)/(1+(b + b1*Genotype_P)*exp(-(c + c1*Genotype_P)*Days));
	OUTPUT OUT = ypred P = yhat;
RUN;
PROC SORT DATA = ypred;
	BY Genotype Days;
RUN;
PROC SGPLOT DATA = ypred;
	SCATTER x = Days y = Leaf_Weight / group = Genotype;
	SERIES    x = Days y = yhat / group = Genotype;
RUN;

* NO a1 --------------------------------------------------------------;
PROC NLIN DATA = soybean_data2;
	PARMS a = 19,  b = 1035, b1 = -490 c = 0.125 c1 = -0.01;
	MODEL Leaf_Weight = (a)/(1+(b + b1*Genotype_P)*exp(-(c + c1*Genotype_P)*Days));
	OUTPUT OUT = ypred P = yhat;
RUN;
PROC SORT DATA = ypred;
	BY Genotype Days;
RUN;
PROC SGPLOT DATA = ypred;
	SCATTER x = Days y = Leaf_Weight / group = Genotype;
	SERIES    x = Days y = yhat / group = Genotype;
RUN;


*CLEARS SAS LOG AND RESULTS FOR CLEANER WORKING ENVIRONMENT;
dm "log; clear; odsresults; clear;";

* RANDOM EFFECTS;
TITLE "Random Effects";   
PROC NLMIXED DATA = soybean_data2;
    PARMS a = 16 ap = 4.78 b = 700 c = 0.125, s2ai = 10, s2bi = 30, s2ci = 0.000001, s2 = 3;
    pred = (a+ap*Genotype_P+ai)/(1+(b+bi)*exp(-(c+ci)*Days));
    MODEL Leaf_Weight~ normal(pred,s2);
    RANDOM ai bi ci ~ normal([0,0,0],[s2ai,0,s2bi,0,0,s2ci]) SUBJECT = Plant_ID;
RUN;

TITLE "Random Effects: Remove ci";   
PROC NLMIXED DATA = soybean_data2;
    PARMS a = 16 ap = 4.78 b = 700 c = 0.125, s2ai = 10, s2bi = 30, s2 = 3;
    pred = (a+ap*Genotype_P+ai)/(1+(b+bi)*exp(-c*Days));
    MODEL Leaf_Weight~ normal(pred,s2);
    RANDOM ai bi ~ normal([0,0],[s2ai,0,s2bi]) SUBJECT = Plant_ID;
RUN;

TITLE "Random Effects: Remove ai";   
PROC NLMIXED DATA = soybean_data2;
    PARMS a = 16 ap = 4.78 b = 90 c = 0.125, s2bi = 30, s2 = 3;
    pred = (a+ap*Genotype_P)/(1+(b+bi)*exp(-c*Days));
    MODEL Leaf_Weight~ normal(pred,s2);
    RANDOM bi ~ normal([0],[s2bi]) SUBJECT = Plant_ID;
RUN;

* RANDOM EFFECTS;
TITLE "Random Effects: ai only";   
ODS SELECT ParameterEstimates;
PROC NLMIXED DATA = soybean_data2;
    PARMS a = 16 ap = 4.78 b = 700 c = 0.125, s2ai = 1, s2 = 3, psi = 0.88;
    pred = (a+ap*Genotype_P+ai)/(1+b*exp(-c*Days));
    MODEL Leaf_Weight~ normal(pred,(pred**(2*psi))*s2);
    RANDOM ai  ~ normal(0,s2ai) SUBJECT = Plant_ID;
RUN;
