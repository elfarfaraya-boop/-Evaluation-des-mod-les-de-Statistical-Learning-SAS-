proc iml;
call randseed(12345);

N = 100;
N_Sim = 1000;
sigma = 0.5;

X_eval = do(-0.9, 0.9, 0.1)`;
Y_true_eval = 2*(X_eval##2) + 1;

Pred_M0 = j(nrow(X_eval), N_Sim, 0);
Pred_M2 = j(nrow(X_eval), N_Sim, 0);
Pred_M8 = j(nrow(X_eval), N_Sim, 0);

do rep = 1 to N_Sim;

    X = j(N, 1);            
    call randgen(X, "Uniform");
    X = X*2 - 1;
   
    Eps = j(N, 1);
    call randgen(Eps, "Normal");
    Eps = Eps * sigma;
   
    Y = 2*(X##2) + 1 + Eps;      

    X_d0 = j(N, 1, 1);

    X_d2 = X_d0 || X || (X##2);

    X_d8 = X_d2;
    do p = 3 to 8; X_d8 = X_d8 || (X##p); end;

    Beta0 = inv(X_d0` * X_d0) * X_d0` * Y;
    Beta2 = inv(X_d2` * X_d2) * X_d2` * Y;
    Beta8 = ginv(X_d8` * X_d8) * X_d8` * Y;

    X_eval_d0 = j(nrow(X_eval), 1, 1);
    X_eval_d2 = X_eval_d0 || X_eval || (X_eval##2);
    X_eval_d8 = X_eval_d2;
    do p = 3 to 8; X_eval_d8 = X_eval_d8 || (X_eval##p); end;

    Pred_M0[,rep] = X_eval_d0 * Beta0;
    Pred_M2[,rep] = X_eval_d2 * Beta2;
    Pred_M8[,rep] = X_eval_d8 * Beta8;
   
end;

Mean_Pred0 = Pred_M0[,:];
Mean_Pred2 = Pred_M2[,:];
Mean_Pred8 = Pred_M8[,:];

BiasSq_0 = (Mean_Pred0 - Y_true_eval)##2;
BiasSq_2 = (Mean_Pred2 - Y_true_eval)##2;
BiasSq_8 = (Mean_Pred8 - Y_true_eval)##2;

Var_0 = j(nrow(X_eval), 1, 0);
Var_2 = j(nrow(X_eval), 1, 0);
Var_8 = j(nrow(X_eval), 1, 0);

do i = 1 to nrow(X_eval);
   Var_0[i] = var(Pred_M0[i,]`);
   Var_2[i] = var(Pred_M2[i,]`);
   Var_8[i] = var(Pred_M8[i,]`);
end;

Resultats = (BiasSq_0[:] || Var_0[:]) //
            (BiasSq_2[:] || Var_2[:]) //
            (BiasSq_8[:] || Var_8[:]);
           
RowNames = {"Degre 0 (Under)", "Degre 2 (Bon)", "Degre 8 (Over)"};
ColNames = {"Biais_Carre", "Variance"};

print "--- COMPROMIS BIAIS VARIANCE (Degrés 0, 2, 8) ---";
print Resultats[rowname=RowNames colname=ColNames];

PlotData = X_eval || Y_true_eval || Pred_M0[,N_Sim] || Pred_M2[,N_Sim] || Pred_M8[,N_Sim];
create PlotDS from PlotData[colname={"X" "True" "Deg0" "Deg2" "Deg8"}];
append from PlotData;
close PlotDS;

quit;

title "Illustration Overfitting (Deg 8) vs Underfitting (Deg 0)";
proc sgplot data=PlotDS;

    series x=X y=True / lineattrs=(color=black thickness=2) legendlabel="Vrai Modèle";
   
    series x=X y=Deg0 / lineattrs=(color=blue pattern=dash) legendlabel="Underfitting (Deg 0)";
   
    series x=X y=Deg2 / lineattrs=(color=green) legendlabel="Bon Modèle (Deg 2)";

    series x=X y=Deg8 / lineattrs=(color=red) legendlabel="Overfitting (Deg 8)";
   
    scatter x=X y=True / markerattrs=(symbol=circlefilled size=5 color=black);
run;
