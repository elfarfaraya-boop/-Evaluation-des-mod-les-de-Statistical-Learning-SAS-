/*Exploratory Data Analysis*/
data WORK.DIABETES;
    infile "/home/u64142266/Modelisation Stoch/Analyse Empirique/Diabetes.txt" 
           dlm='09'x 
           firstobs=3 
           missover;

    input AGE SEX BMI BP S1 S2 S3 S4 S5 S6 Y_char :$20.;
    
    Y = input(compress(Y_char, '\;'), best.);

    drop Y_char;
run;


title "1. Selection STEPWISE avec critere AIC";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=stepwise(select=AIC choose=validate) 
          stats=all; 
run;

title "2. Selection STEPWISE avec critere SBC";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=stepwise(select=SBC choose=validate) 
          stats=all; 
run;


title "3. Regression LASSO";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=lasso(stop=CV choose=CV) 
          cvmethod=random(10); 
run;

title "4. Regression ELASTIC NET";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=elasticnet(stop=CV choose=CV) 
          cvmethod=random(10);
run;

title;

/*CORRELATION MATRIX*/
data DIABETES;
    infile "/home/u64142266/Modelisation Stoch/Analyse Empirique/Diabetes.txt" dlm='09'x firstobs=3;
    input AGE SEX BMI BP S1 S2 S3 S4 S5 S6 Y;
run;


ods graphics on;

proc corr data=DIABETES; 
    var AGE -- S6;
run;


proc iml;
    use DIABETES;
    read all var {AGE SEX BMI BP S1 S2 S3 S4 S5 S6} into MatX;
    close DIABETES;

    R = corr(MatX);
    
    noms = {AGE SEX BMI BP S1 S2 S3 S4 S5 S6};
    print R[colname=noms rowname=noms format=5.2];
quit;

/*BOXPLOT*/

data DIABETES;
    infile "/home/u64142266/Modelisation Stoch/Analyse Empirique/Diabetes.txt" dlm='09'x firstobs=3;
    input AGE SEX BMI BP S1 S2 S3 S4 S5 S6 Y;
run;

proc stdize data=DIABETES out=DIABETES_STD;
    var AGE -- S6;
run;

data DIABETES_PLOT;
    set DIABETES_STD;
    array vars[*] AGE -- S6;
    do i = 1 to dim(vars);
        Variable = vname(vars[i]);
        Value = vars[i];
        output;
    end;
    keep Variable Value;
run;

title "Variable Distribution (Standardized)";
proc sgplot data=DIABETES_PLOT;
    vbox Value / category=Variable group=Variable;
    xaxis display=(nolabel);
    yaxis label="Standardized Value";
run;
title;

/*Simulations*/
data WORK.DIABETES;
    infile "/home/u64142266/Modelisation Stoch/Analyse Empirique/Diabetes.txt" 
           dlm='09'x 
           firstobs=3 
           missover;

    input AGE SEX BMI BP S1 S2 S3 S4 S5 S6 Y_char :$20.;
    
    Y = input(compress(Y_char, '\;'), best.);

    drop Y_char;
run;


title "1. Selection STEPWISE avec critere AIC";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=stepwise(select=AIC choose=validate) 
          stats=all; 
run;

title "2. Selection STEPWISE avec critere SBC";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=stepwise(select=SBC choose=validate) 
          stats=all; 
run;


title "3. Regression LASSO";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=lasso(stop=SBC choose=CV) 
          cvmethod=random(10); 
run;

title "4. Regression ELASTIC NET";
proc glmselect data=WORK.DIABETES plots=(CriterionPanel ASEPlot CoefficientPanel) seed=12345;
    partition fraction(validate=0.3);
    
    model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6 / 
          selection=elasticnet(stop=SBC choose=CV) 
          cvmethod=random(10);
run;

title;
