% Neoclassical Stochastic Growth Model
% Example based on Diaz-Gimenez (2001) - Linear Quadratic Approximations: An Introduction 
% (Chapter 2 in 'Computational Methods for the Study of Dynamic Economies', edited by Marimon & Scott)

% This model is also used by Aldrich, Fernandez-Villaverde, Gallant, & Rubio-Ramirez (2011) - "Tapping the supercomputer under your desk: Solving dynamic equilibrium models with graphics processors,"
% But they do use slightly different parameters to those used here.

% This example solves using anisotropic Smolyak-Chebyshev method
% See Judd, Maliar, Maliar, & Valero (2014) - "Smolyak Method for Solving Dynamic Economic Models: Lagrange Interpolation, Anisotropic Grid and Adaptive Domain"
% Except here the anisotropic Smolyak-Chebyshev method is used as part of value function iteration rather than policy function iteration (on the FOCs)

% NOTE TO SELF: CONSIDER USING q=4 FOR TAUCHEN METHOD. THE DYNARE
% SIMULATION OF 3000 PERIODS GOES OUTSIDE MAX AND MIN VALUES OF Z WHEN
% USING q=3 (am doing this; using q=4)


%% Set up
tauchenoptions.parallel=0; % Use single CPU

Params.beta=0.984;
Params.gamma=2;
Params.alpha=0.35;
Params.delta=0.01;
Params.rho=0.95;
Params.sigma_epsilon=0.005;
Params.sigmasq_epsilon=Params.sigma_epsilon^2;

n_k=501;
n_z=25; % A parameter needed for the Tauchen Method


%% Compute the steady state
K_ss=((Params.alpha*Params.beta)/(1-Params.beta*(1-Params.delta)))^(1/(1-Params.alpha));
X_ss= Params.delta*K_ss;
% These are not really needed; we just use them to determine the grid on
% capital. I mainly calculate them to stay true to original article.

%% Now, create the return function
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val, s_val, gamma, alpha, delta) StochasticNeoClassicalGrowthModel_ReturnFn(aprime_val, a_val, s_val, gamma, alpha, delta);
ReturnFnParamNames={'gamma', 'alpha', 'delta'}; %It is important that these are in same order as they appear in 'StochasticNeoClassicalGrowthModel_ReturnFn'

%% Use Dynare to solve model and get an initial guess of value function.
% Dynare solution using second-order perturbation methods.

% Save the parameters so that dynare can use them to solve 2nd perturbation version of model which is used to get initial guess of value function
beta=Params.beta;
gammac=Params.gamma;
alpha=Params.alpha;
delta=Params.delta;
rho=Params.rho;
sigma_epsilon=Params.sigma_epsilon;
save ./DynareInitialValueFns/paramsfordynare.mat beta gammac alpha delta rho sigma_epsilon

% Since the maintainers of Dynare are a bunch of useless bums who clearly
% don't use matlab and ubuntu the 'dynare-matlab' package is out of date
% and can't be installed. Instead you will just have to manually run dynare
% from Octave after having created 'paramsfordynare.mat', and then come
% back here to load up oo_. and get the policy functions from that.
% dynare ./DynareInitialValueFns/StochasticNeoclassicalGrowthModel.mod

% Load oo_ which contains the policy functions
load ./DynareInitialValueFns/StochasticNeoclassicalGrowthModel_results.mat oo_ M_ options_


% Put the policy functions into a more useable form. Namely make them look
% like the form in which Dynare prints them to screen.
% (the code to do this, next ~100 lines, is a lightly edited copy of Wouter Den Haans 'disp_dr.m')
nx =size(oo_.dr.ghx,2);
nu =size(oo_.dr.ghu,2);
if options_.block
    k = M_.nspred;
    k1 = 1:M_.endo_nbr;
else
    k = find(oo_.dr.kstate(:,2) <= M_.maximum_lag+1);
    klag = oo_.dr.kstate(k,[1 2]);
    k1 = oo_.dr.order_var;
end

var_list = M_.endo_names(1:M_.orig_endo_nbr, :);
nvar = size(M_.endo_names,1);

ivar=zeros(nvar,1);
for i=1:nvar
    i_tmp = strmatch(var_list(i,:),M_.endo_names(k1,:),'exact');
    if isempty(i_tmp)
        error (['One of the variable specified does not exist']) ;
    else
        ivar(i) = i_tmp;
    end
end

decision = [];
for i=1:nvar
    x = oo_.dr.ys(k1(ivar(i)));     
    if options_.order > 1
        x = x + oo_.dr.ghs2(ivar(i))/2;
    end
    decision = [decision x];
end

if options_.order > 1
    decisiontemp = [];
    for i=1:nvar
        x = oo_.dr.ghs2(ivar(i))/2;
        decisiontemp = [decisiontemp x];
    end
    decision = [decision;decisiontemp];
end
% ghx
for k=1:nx
    decisiontemp = [];
    for i=1:nvar
        x = oo_.dr.ghx(ivar(i),k);
        decisiontemp = [decisiontemp x];
    end
    decision = [decision;decisiontemp];
end
% ghu
for k=1:nu
    decisiontemp = [];
    for i=1:nvar
        x = oo_.dr.ghu(ivar(i),k);
        decisiontemp = [decisiontemp x];
    end
    decision = [decision;decisiontemp];
end

if options_.order > 1
    % ghxx
    for k = 1:nx
        for j = 1:k
            decisiontemp = [];
            for i=1:nvar
                if k == j
                    x = oo_.dr.ghxx(ivar(i),(k-1)*nx+j)/2;
                else
                    x = oo_.dr.ghxx(ivar(i),(k-1)*nx+j);
                end
                decisiontemp = [decisiontemp x];                
            end
            decision = [decision;decisiontemp];
        end
    end
    % ghuu
    for k = 1:nu
        for j = 1:k
            decisiontemp = [];
            for i=1:nvar
                if k == j
                    x = oo_.dr.ghuu(ivar(i),(k-1)*nu+j)/2;
                else
                    x = oo_.dr.ghuu(ivar(i),(k-1)*nu+j);
                end
                
                decisiontemp = [decisiontemp x];
            end
            decision = [decision;decisiontemp];
        end
    end
    % ghxu
    for k = 1:nx
        for j = 1:nu
            decisiontemp = [];
            for i=1:nvar
                x = oo_.dr.ghxu(ivar(i),(k-1)*nu+j);
                decisiontemp = [decisiontemp x];
            end
            decision = [decision; decisiontemp];
        end
    end
end
% Policy functions from Dynare (2nd order perturbation) are now contained in 'decision'

%% Now we need to figure out a grid. For capital I take 0.5 times the minimum simulated
% values and 1.5 times the maximum simulated values.
kmin=0.5*min(oo_.endo_simul(2,:));
kmax=1.5*max(oo_.endo_simul(2,:));
k_gridtemp=linspace(kmin,kmax,n_k)'; % Used for initial guess of value fn based on dynare (2nd order perturbation)

% Those for the exogeneous shock the max and min grid pts are taken from Tauchen method using q=4
q=4; % A parameter needed for the Tauchen Method
[z_grid_tauchen, pi_z_tauchen]=TauchenMethod(0,Params.sigmasq_epsilon,Params.rho,n_z,q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q), transmatix is (z,zprime)
zmin=z_grid_tauchen(1);
zmax=z_grid_tauchen(end);

%% Create an initial guess for the value function, first just at values on a fine grid. (also, create version of policy fn for graphing)
% This code is ridiculously slow, but only ever going to be run a few times so who cares.
k_ss_dynare=sum(decision(1:2,2));
z_ss_dynare=sum(decision(1:2,3));
N_epsilon=1; % Note that based on the theory of the model we know that solution depends on z, not z(-1) & epsilon, so just use lots of z(-1) points and no epsilon points gives me just as many z points as having epsilon.

V_Dynare=zeros(length(k_gridtemp),length(z_grid_tauchen));
Vdist=inf;
while Vdist>(10^-9)
    Vold=V_Dynare;
    for k_c=1:length(k_gridtemp)
        k_val=k_gridtemp(k_c);
        k_dev_from_ss=k_val-k_ss_dynare; % k in deviation from steady state
        for z_c=1:length(z_grid_tauchen)
            z_val=z_grid_tauchen(z_c);
            %            z_dev_from_ss=z_val-z_ss_dynare; % z in deviation from steady state
            % Draw an epsilon value, then calculate a z(-1) and hence a k_prime.
            V_kz=nan(N_epsilon,1);
            for epsilon_c=1:N_epsilon % THIS LOOP IS COMPLETELY POINTLESS. THEORY TELLS US THAT WHAT MATTERS IS Z, NOT ZMINUS1 AND EPSILON. SO THIS LOOP JUST (CORRECTLY) GIVES SAME VALUES FOR kprime_val REGARDLESS OF EPSILON
                epsilon_val=Params.sigma_epsilon*randn(1,1);
                zminus1_val=(z_val-epsilon_val)/Params.rho;
                zminus1_dev_from_ss=zminus1_val-z_ss_dynare;
                kprime_val=decision(1,2)+decision(3,2)*k_dev_from_ss+decision(4,2)*zminus1_dev_from_ss+decision(5,2)*epsilon_val+...
                    decision(6,2)*k_dev_from_ss^2+decision(7,2)*k_dev_from_ss*zminus1_dev_from_ss+decision(8,2)*zminus1_dev_from_ss^2+...
                    decision(9,2)*epsilon_val^2+decision(10,2)*k_dev_from_ss*epsilon_val+decision(11,2)*zminus1_dev_from_ss*epsilon_val; % decision
                EV_z=Vold*pi_z_tauchen(z_c,:)';
                % Probably need following line too based on my experience with pure discretized VFI
                EV_z(isnan(EV_z))=0; % Are created when values V=-inf are multiplied by transition probabilities of zero
                V_kz(epsilon_c)=ReturnFn(kprime_val, k_val, z_val, Params.gamma, Params.alpha, Params.delta)+beta*interp1(k_gridtemp, EV_z, kprime_val);
            end
            V_Dynare(k_c,z_c)=sum(V_kz)/N_epsilon;
        end
    end
    Vdist=max(reshape(abs(V_Dynare-Vold),[numel(V_Dynare),1]))
end
    
% Now the policy fn
Policy_Dynare=zeros(length(k_gridtemp),length(z_grid_tauchen));
for k_c=1:length(k_gridtemp)
    k_val=k_gridtemp(k_c);
    k_dev_from_ss=k_val-k_ss_dynare; % k in deviation from steady state
    for z_c=1:length(z_grid_tauchen)
        z_val=z_grid_tauchen(z_c);
        Policy_kz=nan(N_epsilon,1);
        %         for epsilon_c=1:N_epsilon
        epsilon_val=0; %Params.sigma_epsilon*randn(1,1);
        zminus1_val=(z_val-epsilon_val)/Params.rho;
        zminus1_dev_from_ss=zminus1_val-z_ss_dynare;
        kprime_val=decision(1,2)+decision(3,2)*k_dev_from_ss+decision(4,2)*zminus1_dev_from_ss+decision(5,2)*epsilon_val+...
            decision(6,2)*k_dev_from_ss^2+decision(7,2)*k_dev_from_ss*zminus1_dev_from_ss+decision(8,2)*zminus1_dev_from_ss^2+...
            decision(9,2)*epsilon_val^2+decision(10,2)*k_dev_from_ss*epsilon_val+decision(11,2)*zminus1_dev_from_ss*epsilon_val; % decision
        %             Policy_kz(epsilon_c)=kprime_val;
        %         end
        Policy_Dynare(k_c,z_c)=kprime_val; %mean(Policy_kz);
    end
end

save ./SavedOutput/StochNeoClassicalGrowth_VfromDynare.mat V_Dynare Policy_Dynare k_gridtemp z_grid_tauchen

% Check if V contains -Inf values as Smolyak Chebyshev cannot handle these.
if min(min(V_Dynare))==-Inf
    fprintf('WARNING: VALUE FN CONTAINS A -INF VALUE. THIS WILL CAUSE PROBLEMS FOR SMOLYAK-CHEBYSHEV APPROXIMATION \n')
end

%% Plot the Dynare policy and value function

[x1,x2] = meshgrid(z_grid_tauchen,k_gridtemp);

figure(4)
surf(x1, x2, V_Dynare)
title('Value Fn: Dynare')

figure(5)
surf(x1,x2,Policy_Dynare)
title('Policy Fn: Dynare')




