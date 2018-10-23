%% Using Smolyak-Chebyshev to fit the Value Function itself. What goes wrong

% Before running this you must run the following which creates V_Dynare and
% Policy_Dynare (and in turn requires that you run a Dynare mod-file)
% StochasticNeoClassicalGrowthModel_SmolyakChebyshev2.m

%% 1. Smolyak anisotropic method for 2 dimensions: set-up
% -----------------------------------------------
d1_orders=[1,2,3,4,5,7,9];
d2_orders=[1,2,3,4,5,7,9];
n1=length(d1_orders);
n2=length(d2_orders);
L2_V_Table=zeros(n1,n2);
LInf_V_Table=zeros(n1,n2);
L2_Policy_Table=zeros(n1,n2);
LInf_Policy_Table=zeros(n1,n2);

for ii=1:n1
    for jj=1:n2
        vector_mus_dimensions = [d1_orders(ii),d2_orders(jj)]; % Introduce the level of approximation in every dimension from 1 to 10; see Section 4 of JMMV (2014)
        [L2_V,LInf_V, L2_Policy, LInf_Policy]=WhatGoesWrong_Table_fn(vector_mus_dimensions);
        L2_V_Table(ii,jj)=L2_V;
        LInf_V_Table(ii,jj)=LInf_V;
        L2_Policy_Table(ii,jj)=L2_Policy;
        LInf_Policy_Table(ii,jj)=LInf_Policy;
    end
end

% Look at some outcomes
LInf_V_Table(6:7,6:7)
LInf_Policy_Table(6:7,6:7)

L2_V_Table(6:7,6:7)
L2_Policy_Table(6:7,6:7)

LInf_V_Table(2,2)
LInf_Policy_Table(2,2)

%Table
FID = fopen('./FiguresForPaper/Table_ValueFn_LInf.tex', 'w');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l|ccccccc} \n \\hline \\hline \n');
fprintf(FID, ' Order \\textbackslash Order & %d & %d & %d & %d & %d & %d & %d \\\\ \n \\hline \n',d2_orders);
for ii=1:length(d1_orders)
    fprintf(FID, '%d & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e \\\\ \n ', d1_orders(ii),LInf_V_Table(ii,:));
end
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Notes: $L_{\\infty}$ norm measure of difference between Value function (for Stochastic Neoclassical Growth Model) and the Smolyak-Chebyshev approximation of given orders. \\\\ \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%Table
FID = fopen('./FiguresForPaper/Table_PolicyFn_LInf.tex', 'w');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}l|ccccccc} \n \\hline \\hline \n');
fprintf(FID, ' Order \\textbackslash Order & %d & %d & %d & %d & %d & %d & %d \\\\ \n \\hline \n',d2_orders);
for ii=1:length(d1_orders)
    fprintf(FID, '%d & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e & %8.2e \\\\ \n ', d1_orders(ii),LInf_Policy_Table(ii,:));
end
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Notes: $L_{\\infty}$ norm measure of difference between optimal policy function (for Stochastic Neoclassical Growth Model) and the Smolyak-Chebyshev approximation of given orders. \\\\ \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

