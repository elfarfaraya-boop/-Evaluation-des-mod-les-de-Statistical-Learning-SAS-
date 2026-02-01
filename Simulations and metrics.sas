/*DGP GAUSSIEN*/
proc datasets lib=work nolist;
    delete DGP_Gaussien Tableau_Recapitulatif ParEst Statistiques Recap_Ligne;
quit;

proc iml;

call randseed(12345);
n  = 200;
p  = 50;
MC = 1000;

Beta = j(p, 1, 0);
Beta[1:6] = {1, -0.35, 0.15, 0.27, 0.57, -0.14};

Mu       = j(1, p, 0);
Varcovar = I(p);

simulated_data = j(n*MC, 2+p, .);
a = 1;

do iteration = 1 to MC;
  X   = RandNormal(n, Mu, Varcovar);
  eps = normal(j(n, 1, 0))*0.1;
  Y   = X * Beta + eps;

  simulated_data[a:a+n-1, 1]     = iteration;
  simulated_data[a:a+n-1, 2]     = Y;
  simulated_data[a:a+n-1, 3:2+p] = X;
  a = a + n;
end;

cname = {"Iteration_ID" "Y"};
do i = 1 to p;
  cname = cname || cats("X", i);
end;

create DGP_Gaussien from simulated_data[colname=cname];
append from simulated_data;
close DGP_Gaussien;

quit;

/*Multicollinearity Simulation*/
proc iml;
call randseed(12345);

N  = 200;
P  = 50;
MC = 1000;

Beta = j(P, 1, 0);
Beta[1:6] = {1, -0.35, 0.15, 0.27, 0.57, -0.14};

Mu = j(1, P, 0);

r = {1, 0.8, 0.75, 0.7, 0.65, 0.6};   
V1 = toeplitz(r);                
Varcovar = I(P);
Varcovar[1:6, 1:6] = V1;

simulated_data = j(N*MC, 2+P, .);
a = 1;

do iteration = 1 to MC;
	X = RandNormal(N, Mu, Varcovar);

    eps = normal(j(N,1,0)) * 0.1;
    Y   = X * Beta + eps;

    simulated_data[a:a+N-1, 1]     = iteration;
    simulated_data[a:a+N-1, 2]     = Y;
    simulated_data[a:a+N-1, 3:2+P] = X;

    a = a + N;
end;

cname = {"Iteration_ID" "Y"};
do i = 1 to P;
    cname = cname || cats("X", i);
end;

create DGP_Multicolinearite from simulated_data[colname=cname];
append from simulated_data;
close DGP_Multicolinearite;

quit;

/*Outliers simulation*/
proc iml;
start Outliers(X_in);
	X_out = X_in;
	n_obs = nrow(X_in);
    p_var = ncol(X_in);


	do i = 1 to n_obs ;
		do j = 1 to p_var;
			u = uniform(0); 
			if u > 0.90 then do;
				X_out[i,j] = 4 + normal(0);	 
			end;
		end;
	end;
	return(X_out);
finish;

store module=(Outliers);
load module=( Outliers);
call randseed(12345);
    
N = 200;
P = 50; 
MC = 1000;

Mu = j(1, p, 0);      
Varcovar = I(p);    

Beta = j(P, 1, 0); 
Beta[1:6] = {1, -0.35, 0.15, 0.27, 0.57, -0.14};

simulated_data = j(N*MC, 2+P, .);
a = 1;

do iteration = 1 to MC;

    X_0     = RandNormal(N, Mu, Varcovar);
    X_final = Outliers(X_0);

    eps = normal(j(N, 1, 0)) * 0.1;
    Y   = X_final * Beta + eps;

    simulated_data[a:a+N-1, 1]     = iteration;
    simulated_data[a:a+N-1, 2]     = Y;
    simulated_data[a:a+N-1, 3:2+P] = X_final;

    a = a + N;
end;

cname = {"Iteration_ID" "Y"};
do i = 1 to P;
    cname = cname || cats("X", i);
end;

create DGP_Outliers from simulated_data[colname=cname];
append from simulated_data;
close DGP_Outliers;

quit;

/*Structural break simulation*/

proc iml;
call randseed(12345);
n  = 200;
P  = 50;
MC = 1000;

Beta1 = j(P, 1, 0);
Beta1[1:6] = {1, -0.35, 0.15, 0.27, 0.57, -0.14};          

Beta2 = j(P, 1, 0);
Beta2[1:7] = {-1, 0.55, 0.6, -0.35, 0.15, 0.45, 0.1};   

Mu       = j(1, P, 0);
Varcovar = I(P);

simulated_data = j(n*MC, 2+P, .);
a = 1;

do iteration = 1 to MC;

    t_break = rand("Integer", floor(0.1*n), floor(0.7*n));

    X   = RandNormal(n, Mu, Varcovar);
    eps = normal(j(n, 1, 0)) * 0.1;
    Y   = j(n, 1, .);

    Y[1:t_break]           = X[1:t_break, ] * Beta1 + eps[1:t_break];
    Y[(t_break+1):n]       = X[(t_break+1):n, ] * Beta2 + eps[(t_break+1):n];

    simulated_data[a:a+n-1, 1]     = iteration;
    simulated_data[a:a+n-1, 2]     = Y;
    simulated_data[a:a+n-1, 3:2+P] = X;

    a = a + n;
end;

cname = {"Iteration_ID" "Y"};
do i = 1 to P;
    cname = cname || cats("X", i);
end;

create DGP_Rupture from simulated_data[colname=cname];
append from simulated_data;
close DGP_Rupture;
quit;

/*Mix break/outliers simulation*/


proc iml;
    
    load module=(Outliers);

    call randseed(12345);

    
    n  = 200;
    P  = 50;
    MC = 1000;

    Mu       = j(1, P, 0);
    Varcovar = I(P);

    
    Beta1 = j(P, 1, 0);
    Beta1[1:6] = {1, -0.35, 0.15, 0.27, 0.57, -0.14};          

    Beta2 = j(P, 1, 0);
    Beta2[1:7] = {-1, 0.55, 0.6, -0.35, 0.15, 0.45, 0.1};     

    simulated_data = j(n*MC, 2+P, .);
    a = 1;

    do iteration = 1 to MC;

        t_break = rand("Integer", floor(0.1*n), floor(0.7*n));

        X0  = RandNormal(n, Mu, Varcovar);
        X   = Outliers(X0);

        eps = normal(j(n, 1, 0)) * 0.1;
        Y   = j(n, 1, .);

        Y[1:t_break]     = X[1:t_break, ] * Beta1 + eps[1:t_break];
        Y[(t_break+1):n] = X[(t_break+1):n, ] * Beta2 + eps[(t_break+1):n];

        simulated_data[a:a+n-1, 1]     = iteration;
        simulated_data[a:a+n-1, 2]     = Y;
        simulated_data[a:a+n-1, 3:2+P] = X;

        a = a + n;
    end;

    cname = {"Iteration_ID" "Y"};
    do i = 1 to P;
        cname = cname || cats("X", i);
    end;

    create DGP_BreakOutliers from simulated_data[colname=cname];
    append from simulated_data;
    close DGP_BreakOutliers;

quit;

/*Mix multicollinearity/outliers simulation*/

proc iml;
start ImanConover(X, C);
n_obs = nrow(X);
p_var = ncol(X);
S = J(n_obs, p_var);

do i = 1 to ncol(X);
ranks = ranktie(X[,i], "mean");
S[,i] = quantile("Normal", ranks/(n_obs+1));
end;

CS = corr(S); 
Q = root(CS);
P = root(C); 

T = solve(Q,P); 
Y = S*T; 

W = X;
do i = 1 to ncol(Y);
r_vec = rank(Y[,i]); 
tmp = W[,i];
call sort(tmp); 
W[,i] = tmp[r_vec]; 
end;

return( W );
finish;

store module=(ImanConover);

load module=(ImanConover Outliers);

call randseed(12345);

N  = 200;
P  = 50;
MC = 1000;

Mu       = j(1, P, 0);
Varcovar = I(P);

Beta = j(P, 1, 0);
Beta[1:6] = { 1, -0.35, 0.15, 0.27, 0.57, -0.14 };

rho = {1, 0.8, 0.75, 0.7, 0.65, 0.6};
C = toeplitz(rho);

simulated_data = j(N*MC, 2+P, .);
a = 1;

do iteration = 1 to MC;

	X_indep = RandNormal(N, Mu, Varcovar);

	X_outl = Outliers(X_indep);

	X1 = ImanConover(X_outl[, 1:6], C);
    X_final = X1 || X_outl[, 7:P];

	eps = normal(j(N, 1, 0)) * 0.1;
	Y   = X_final * Beta + eps;

	simulated_data[a:a+N-1, 1]     = iteration;
	simulated_data[a:a+N-1, 2]     = Y;
	simulated_data[a:a+N-1, 3:2+P] = X_final;

	a = a + N;
	end;

cname = {"Iteration_ID" "Y"} || ("X" + strip(char(1:P)));
create DGP_CorrOutliers from simulated_data[colname=cname];
append from simulated_data;
close DGP_CorrOutliers;

quit;

/*Mix multicollinearity/break simulation*/
proc iml;
start ImanConover(X, C);
n_obs = nrow(X);
p_var = ncol(X);
S = J(n_obs, p_var);

do i = 1 to ncol(X);
ranks = ranktie(X[,i], "mean");
S[,i] = quantile("Normal", ranks/(n_obs+1));
end;

CS = corr(S); 
Q = root(CS); 
P = root(C);
T = solve(Q,P); 
Y = S*T; 

W = X;
do i = 1 to ncol(Y);
r_vec = rank(Y[,i]);
tmp = W[,i];
call sort(tmp); 
W[,i] = tmp[r_vec]; 
end;

return( W );
finish;

load module=(ImanConover);

call randseed(12345);

N  = 200;
P  = 50;
MC = 1000;

Mu       = j(1, P, 0);
Varcovar = I(P);

Beta1 = j(P, 1, 0);
Beta1[1:6] = {1, -0.35, 0.15, 0.27, 0.57, -0.14};        

Beta2 = j(P, 1, 0);
Beta2[1:7] = {-1, 0.55, 0.6, -0.35, 0.15, 0.45, 0.1}; 
rho = {1, 0.8, 0.75, 0.7, 0.65, 0.6};
C = toeplitz(rho);

simulated_data = j(N*MC, 2+P, .);
a = 1;

do iteration = 1 to MC;
	
	t_break = rand("Integer", floor(0.1*n), floor(0.7*n));

	X_indep = RandNormal(N, Mu, Varcovar);
	X_corr = ImanConover(X_indep[, 1:6], C);
    X_final = X_corr || X_indep[, 7:P];

    eps = normal(j(N, 1, 0)) * 0.1;
	
	Y   = j(n, 1, .);

    Y[1:t_break]     = X_final[1:t_break, ] * Beta1 + eps[1:t_break];
	Y[(t_break+1):n] = X_final[(t_break+1):n, ] * Beta2 + eps[(t_break+1):n];

	simulated_data[a:a+N-1, 1]     = iteration;
	simulated_data[a:a+N-1, 2]     = Y;
	simulated_data[a:a+N-1, 3:2+P] = X_final;

	a = a + N;
	end;

cname = {"Iteration_ID" "Y"} || ("X" + strip(char(1:P)));
create DGP_CorrBreak from simulated_data[colname=cname];
append from simulated_data;
close DGP_CorrBreak;

quit;

/*Combined scenario*/
proc iml;
start ImanConover(X, C);
n_obs = nrow(X);
p_var = ncol(X);
S = J(n_obs, p_var);

do i = 1 to ncol(X);
ranks = ranktie(X[,i], "mean");
S[,i] = quantile("Normal", ranks/(n_obs+1));
end;

CS = corr(S); 
Q = root(CS); 
P = root(C); 

T = solve(Q,P); 
Y = S*T; 

W = X;
do i = 1 to ncol(Y);
r_vec = rank(Y[,i]); 
tmp = W[,i];
call sort(tmp); 
W[,i] = tmp[r_vec]; 
end;

return( W );
finish;

store module=(ImanConover);

load module=(ImanConover Outliers);

call randseed(12345);

N  = 200;
P  = 50;
MC = 1000;

Mu       = j(1, P, 0);
Varcovar = I(P);

Beta1 = j(P, 1, 0);
Beta1[1:6] = {1, -0.35, 0.15, 0.27, 0.57, -0.14};          

Beta2 = j(P, 1, 0);
Beta2[1:7] = {-1, 0.55, 0.6, -0.35, 0.15, 0.45, 0.1}; 

rho = {1, 0.8, 0.75, 0.7, 0.65, 0.6};
C = toeplitz(rho);

simulated_data = j(N*MC, 2+P, .);
a = 1;

do iteration = 1 to MC;

	t_break = rand("Integer", floor(0.1*n), floor(0.7*n));

	X_indep = RandNormal(N, Mu, Varcovar);
	X_outl = Outliers(X_indep);  
	X1 = ImanConover(X_outl[, 1:6], C);
    X_final = X1 || X_outl[, 7:P];

	eps = normal(j(N, 1, 0)) * 0.1;
	
	Y   = j(n, 1, .);

	Y[1:t_break]     = X_final[1:t_break, ] * Beta1 + eps[1:t_break];
	Y[(t_break+1):n] = X_final[(t_break+1):n, ] * Beta2 + eps[(t_break+1):n];

	simulated_data[a:a+N-1, 1]     = iteration;
	simulated_data[a:a+N-1, 2]     = Y;
	simulated_data[a:a+N-1, 3:2+P] = X_final;

	a = a + N;
	end;

cname = {"Iteration_ID" "Y"} || ("X" + strip(char(1:P)));
create DGP_CorrOutliersBreak from simulated_data[colname=cname];
append from simulated_data;
close DGP_CorrOutliersBreak;

quit;

/*Metric Macro for the Gaussian, Multicollinearity, and Outliers DGPs*/

%macro Metrique(data_in, methode, choose, stop, sle_val=0.15);

  ods exclude all;
  ods output ParameterEstimates = ParEst;

  proc glmselect data=&data_in. plots=none;
      by Iteration_ID;
      model Y = X1-X50 / selection=&methode. (choose=&choose. stop=&stop. sle=&sle_val.);
  run;

  ods exclude none;

  proc iml;
      vrai = {"X1" "X2" "X3" "X4" "X5" "X6"};

      use ParEst;
          read all var {"Iteration_ID"} into IDs;
          read all var {"Parameter"} into Choisi;
      close;

      uniqueIDs = unique(IDs);
      nSim = ncol(uniqueIDs);

      results = j(nSim, 3, .);  

      do i = 1 to nSim;
          currentID = uniqueIDs[i];

          idx = loc(IDs = currentID & upcase(Choisi) ^= "INTERCEPT");

          if ncol(idx) > 0 then do;
              vars_selectionnees   = upcase(Choisi[idx]`);
              Faux_Positifs_list   = setdif(vars_selectionnees, vrai);
              Faux_Negatifs_list   = setdif(vrai, vars_selectionnees);
          end;
          else do;
              vars_selectionnees   = { };
              Faux_Positifs_list   = { };
              Faux_Negatifs_list   = vrai;
          end;

          results[i, 1] = currentID;
          results[i, 2] = ncol(Faux_Positifs_list);
          results[i, 3] = ncol(Faux_Negatifs_list);
      end;

      create Statistiques from results[colname={"Iteration_ID" "Faux_Positifs" "Faux_Negatifs"}];
      append from results;
      close Statistiques;
  quit;

  proc sql;
      create table Recap_Ligne as
      select
          "&methode" as Methode length=20,
          "&choose"  as Choose  length=20,
          case
              when "&stop" = "SL" then "SL=" || put(&sle_val, 5.2)
              else "&stop"
          end as Stop length=20,

          sum(case when Faux_Positifs=0 and Faux_Negatifs=0 then 1 else 0 end) / count(*) as Perfectfit format=percent8.2,
          sum(case when Faux_Positifs>0 and Faux_Negatifs=0 then 1 else 0 end) / count(*) as Overfit    format=percent8.2,
          sum(case when Faux_Positifs=0 and Faux_Negatifs>0 then 1 else 0 end) / count(*) as Underfit   format=percent8.2,
          sum(case when Faux_Positifs>0 and Faux_Negatifs>0 then 1 else 0 end) / count(*) as Wrong      format=percent8.2
      from Statistiques;
  quit;

  proc append base=Tableau_Recapitulatif data=Recap_Ligne force;
  run;

%mend;

/*Metric Macro for the structural break DGP*/
%macro Metrique(data_in, methode, choose, stop, sle_val=0.15);

  ods exclude all;
  ods output ParameterEstimates = ParEst;

  proc glmselect data=&data_in. plots=none;
      by Iteration_ID;
      model Y = X1-X50 / selection=&methode. (choose=&choose. stop=&stop. sle=&sle_val.);
  run;

  ods exclude none;

  proc iml;
      vrai = {"X1" "X2" "X3" "X4" "X5" "X6" "X7"};

      use ParEst;
          read all var {"Iteration_ID"} into IDs;
          read all var {"Parameter"} into Choisi;
      close;

      uniqueIDs = unique(IDs);
      nSim = ncol(uniqueIDs);

      results = j(nSim, 3, .);  

      do i = 1 to nSim;
          currentID = uniqueIDs[i];

          idx = loc(IDs = currentID & upcase(Choisi) ^= "INTERCEPT");

          if ncol(idx) > 0 then do;
              vars_selectionnees   = upcase(Choisi[idx]`);
              Faux_Positifs_list   = setdif(vars_selectionnees, vrai);
              Faux_Negatifs_list   = setdif(vrai, vars_selectionnees);
          end;
          else do;
              vars_selectionnees   = { };
              Faux_Positifs_list   = { };
              Faux_Negatifs_list   = vrai;
          end;

          results[i, 1] = currentID;
          results[i, 2] = ncol(Faux_Positifs_list);
          results[i, 3] = ncol(Faux_Negatifs_list);
      end;

      create Statistiques from results[colname={"Iteration_ID" "Faux_Positifs" "Faux_Negatifs"}];
      append from results;
      close Statistiques;
  quit;

  proc sql;
      create table Recap_Ligne as
      select
          "&methode" as Methode length=20,
          "&choose"  as Choose  length=20,
          case
              when "&stop" = "SL" then "SL=" || put(&sle_val, 5.2)
              else "&stop"
          end as Stop length=20,

          sum(case when Faux_Positifs=0 and Faux_Negatifs=0 then 1 else 0 end) / count(*) as Perfectfit format=percent8.2,
          sum(case when Faux_Positifs>0 and Faux_Negatifs=0 then 1 else 0 end) / count(*) as Overfit    format=percent8.2,
          sum(case when Faux_Positifs=0 and Faux_Negatifs>0 then 1 else 0 end) / count(*) as Underfit   format=percent8.2,
          sum(case when Faux_Positifs>0 and Faux_Negatifs>0 then 1 else 0 end) / count(*) as Wrong      format=percent8.2
      from Statistiques;
  quit;

  proc append base=Tableau_Recapitulatif data=Recap_Ligne force;
  run;

%mend;
