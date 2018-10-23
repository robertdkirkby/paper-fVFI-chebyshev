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

% A few lines needed for running on the Server
addpath(genpath('./MatlabToolkits/'))

SkipStuff=1

Javier=0;   %If you set this parameter to 0 then the parameters will all be set to those used by Aldrich, Fernandez-Villaverde, Gallant, & Rubio-Ramirez (2011)

%% Set up
tauchenoptions.parallel=0; % Use single CPU
vfoptions.parallel=0; % Use single CPU

%Discounting rate
beta = 0.96;

%Give the parameter values (Params will be a 'structure' containing all the parameter values)
Params.alpha = 0.33;
Params.gamma=1; %gamma=1 is log-utility
Params.rho = 0.95;
Params.delta = 0.10;
Params.sigmasq_epsilon=0.09;

if Javier==0
    n_z=4;
    Params.beta=0.984;
    Params.gamma=2;
    Params.alpha=0.35;
    Params.delta=0.01;
    Params.rho=0.95;
    Params.sigma_epsilon=0.005;
    Params.sigmasq_epsilon=Params.sigma_epsilon^2;
    vfoptions.tolerance=(1-beta)*10^(-8); % Toolkit default is 10^(-9)
    vfoptions.howards=20; % Toolkit default is 80
end

%% Compute the steady state
K_ss=((Params.alpha*Params.beta)/(1-Params.beta*(1-Params.delta)))^(1/(1-Params.alpha));
X_ss= Params.delta*K_ss;
%These are not really needed; we just use them to determine the grid on
%capital. I mainly calculate them to stay true to original article.

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

% 1001 for k, 25 for z, 100 for epsilon.

%% Now we need to figure out a grid. For capital I take 0.5 times the minimum simulated
% values and 1.5 times the maximum simulated values.
kmin=0.5*min(oo_.endo_simul(2,:));
kmax=1.5*max(oo_.endo_simul(2,:));
k_gridtemp=linspace(kmin,kmax,501)'; % Used for initial guess of value fn based on dynare (2nd order perturbation)

% Those for the exogeneous shock the max and min grid pts are taken from Tauchen method using q=4
n_z=5; % A parameter needed for the Tauchen Method
q=4; % A parameter needed for the Tauchen Method
[z_grid_tauchen, pi_z_tauchen]=TauchenMethod(0,Params.sigmasq_epsilon,Params.rho,n_z,q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q), transmatix is (z,zprime)
zmin=z_grid_tauchen(1);
zmax=z_grid_tauchen(end);

save ./GPU_Solution/gridendpoints.mat kmin kmax k_gridtemp zmin zmax z_grid_tauchen

%% Create an initial guess for the value function, first just at values on a fine grid. (also, create version of policy fn for graphing)
k_ss_dynare=sum(decision(1:2,2));
z_ss_dynare=sum(decision(1:2,3));
N_epsilon=20; %100

if SkipStuff==0
    V_Dynare=zeros(length(k_gridtemp),length(z_grid_tauchen));
    Vdist=inf;
    while Vdist>(10^-7)
        Vold=V_Dynare;
        for k_c=1:length(k_gridtemp)
            k_val=k_gridtemp(k_c);
            k_dev_from_ss=k_val-k_ss_dynare; % k in deviation from steady state
            parfor z_c=1:length(z_grid_tauchen)
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
        parfor z_c=1:length(z_grid_tauchen)
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
    
    save ./DynareInitialValueFns/StochNeoClassicalGrowth_VfromDynare.mat V_Dynare Policy_Dynare
else
    load ./DynareInitialValueFns/StochNeoClassicalGrowth_VfromDynare.mat V_Dynare Policy_Dynare
end

% Check if V contains -Inf values as Smolyak Chebyshev cannot handle these.
if min(min(V_Dynare))==-Inf
    fprintf('WARNING: VALUE FN CONTAINS A -INF VALUE. THIS WILL CAUSE PROBLEMS FOR SMOLYAK-CHEBYSHEV APPROXIMATION \n')
end

%% Plot the Dynare policy and value function
load ./DynareInitialValueFns/StochNeoClassicalGrowth_VfromDynare.mat V_Dynare

[x1,x2] = meshgrid(z_grid_tauchen,k_gridtemp);

figure(1)
surf(x1, x2, V_Dynare)
title('Value Fn: Dynare')

figure(3)
surf(x1,x2,Policy_Dynare)
title('Policy Fn: Dynare')


%% Create Smolyak-Chebyshev grids (grids are defined as a column vectors)

% Using Smolyak-Chebyshev method to approximate the value function and policy function.
vfoptions.solnmethod='smolyak_chebyshev';
% Set 'level of approximation' for each dimension of the value function (must be given for (d,a,z)).
vfoptions.aniso_levelofapprox=[5,3]; % If you don't set this it is assumed you want '5' for each and every dimension.

k_grid=[kmin; kmax]; % With Smolyak-Chebyshev just give min and max values for the domain
n_k=2;
% put in vfi toolkit internal notation
n_a=n_k;
a_grid=k_grid;

z_grid=[zmin; zmax]; % Essentially: =[z_grid_tauchen(1); z_grid_tauchen(end)]; 
n_z=2;
% To take advantage of the sparse grid on the (continuous) z variables we
% need to set pi_z up a transition matrix between the sparse grid points.
% Creating pi_z for smolyak sparse grid on the continuous z variables is a
% complicated step. We first need to figure out what size it is.
if n_z(1)==0
    numdimensions=length(n_a);
else
    numdimensions=length(n_a)+length(n_z);
end
% Use JMMV2014 codes to fit the Smolyak-Chebyshev approximation
if isfield(vfoptions,'aniso_levelofapprox')==1
    % Introduce the level of approximation in every dimension from 1 to 10; see Section 4 of JMMV (2014)
    vector_mus_dimensions = vfoptions.aniso_levelofapprox; 
else
    vector_mus_dimensions = 5*numdimensions; % Default value
end
mu_max  = max(vector_mus_dimensions);
Smolyak_elem_iso = Smolyak_Elem_Isotrop(numdimensions,mu_max);
Smol_elem_ani = Smolyak_Elem_Anisotrop(Smolyak_elem_iso,vector_mus_dimensions);
Smol_grid_ani = Smolyak_Grid(numdimensions,mu_max,Smol_elem_ani);


%% Create transition matrix for taking expectations.
% We want to calculate the number of unique z grid points. This will tell
% us the required size of pi_z
% First, get all the z values from the anisotropic Smolyak grid on (a,z): Smol_grid_ani
Smol_grid_ani_z = Smol_grid_ani(:,length(n_a)+1:end);
% Find all the unique rows (code deals with general case, here it is
% overkill in sense that there is only one column).
Smol_grid_ani_z_unique=Smol_grid_ani_z(1,:);
jj=1;
for ii=2:size(Smol_grid_ani_z,1)
    temp=Smol_grid_ani_z(ii,:);
    anewone=1;
    for kk=1:jj
        if sum(temp==Smol_grid_ani_z_unique(kk,:))==length(n_z)
            anewone=0;
        end
    end
    if anewone==1
        Smol_grid_ani_z_unique=[Smol_grid_ani_z_unique;Smol_grid_ani_z(ii,:)];
        jj=jj+1;
    end
end
% jj is now the number of unique rows
% So we know that size of pi_z should be
pi_z_unique=zeros(jj,jj); % jj is equal to size(Smol_grid_ani_z_unique,1)
% Now we just need to give it some values that make it a transition matrix (rows sum to 1)
% The present case is easy enough since we want it to be approximating an AR(1).
% So just do a Tauchen style approach to calculating the transition probabilities.
% First though, our z values are currently on the [-1,1] hypercube. Move
% them to the relevant [a,b] z_grid.
Smol_grid_ani_z_unique_ab=(Smol_grid_ani_z_unique*(z_grid(2)-z_grid(1))+z_grid(1)+z_grid(2))/2;
rho=Params.rho;
mew=0;
sigma=sqrt(Params.sigmasq_epsilon);
Smol_grid_ani_z_unique_ab=sort(Smol_grid_ani_z_unique_ab);
for ii=1:length(Smol_grid_ani_z_unique_ab)
    omega1=Smol_grid_ani_z_unique_ab(2)-Smol_grid_ani_z_unique_ab(1);
    pi_z_unique(ii,1)=normcdf(Smol_grid_ani_z_unique_ab(1)+omega1/2-rho*Smol_grid_ani_z_unique_ab(ii),mew,sigma);
    for jj=2:length(Smol_grid_ani_z_unique_ab)-1
        omega1=Smol_grid_ani_z_unique_ab(jj+1)-Smol_grid_ani_z_unique_ab(jj);
        omega2=Smol_grid_ani_z_unique_ab(jj)-Smol_grid_ani_z_unique_ab(jj-1);
        pi_z_unique(ii,jj)=normcdf(Smol_grid_ani_z_unique_ab(jj)+omega1/2-rho*Smol_grid_ani_z_unique_ab(ii),mew,sigma)-normcdf(Smol_grid_ani_z_unique_ab(jj)-omega2/2-rho*Smol_grid_ani_z_unique_ab(ii),mew,sigma);
    end
    omega2=Smol_grid_ani_z_unique_ab(end)-Smol_grid_ani_z_unique_ab(end-1);
    pi_z_unique(ii,end)=1-normcdf(Smol_grid_ani_z_unique_ab(end)-omega2/2-rho*Smol_grid_ani_z_unique_ab(ii),mew,sigma);
end
% % (given that Tauchen method is normally implemented with matrices not for loops this could certainly be improved)
% % For more general case, maybe a simulation combined with nearest
% % neighbour to assign each simulated observation with a point, then divide by number 
% % by number of simulations (do for each starting point). This is not going
% % to be the fastest way, but will work.
% %
% % Now that we have the pi_z_unique, we just create the full pi_Smol_grid_ani which
% % gives the transition probabilities between the various rows of Smol_grid_ani and
% % populate it using the pi_z_unique values.
% pi_Smol_grid_ani=zeros(size(Smol_grid_ani,1),size(Smol_grid_ani,1));
% for ii=1:size(Smol_grid_ani,1) % This for-loop would be much simpler if I had stored the 'unique' indexes back when I created 'unique' from Smol_grid_ani_z
%     for ii2=1:size(Smol_grid_ani_z_unique,1)
%         if Smol_grid_ani_z(ii,:)==Smol_grid_ani_z_unique(ii2,:)
%             ii_c=ii2;
%         end
%     end
%     for jj=1:size(Smol_grid_ani,1)
%         for jj2=1:size(Smol_grid_ani_z_unique,1)
%             if Smol_grid_ani_z(jj,:)==Smol_grid_ani_z_unique(jj2,:)
%                 jj_c=jj2;
%             end
%         end
%         % Smol_grid_ani has both 'a' and 'z'. We have found the 'z' and
%         % 'zprime' and we have the probability of moving between stored in
%         % pi_z_unique. If both of the 'a' are the same then we want to
%         % store this value in pi_Smol_grid_ani.
%         if prod(Smol_grid_ani(ii,1:length(n_a))==Smol_grid_ani(jj,1:length(n_a)))==1 % If all of the 'a' elements are equal
%             pi_Smol_grid_ani(ii,jj)=pi_z_unique(ii_c,jj_c);
%         end
%     end
% end % This implementation is highly inefficient. No need for them to be nested in this way. (will only run once so can't be bothered fixing)
% for ii=1:size(Smol_grid_ani,1)
%     pi_Smol_grid_ani(ii,:)=pi_Smol_grid_ani(ii,:)/sum(pi_Smol_grid_ani(ii,:)); % Normalize rows so that each row sums to one.
% end % This loop is unnecessary, could be done as matrix operation
% %
% % Phew! Done. Now have pi_z for a sparse Smolyak-Chebyshev grid!!!
% % surf(cumsum(pi_Smol_grid_ani,2))


%% Create transition matrix for taking expectations: alternative method.
% mew=0; sigma=sigma_epsilon;
% % we just create the full pi_Smol_grid_ani which
% % gives the transition probabilities between the various rows of Smol_grid_ani and
% % populate it.
% % First step is to create a binary matrix that just matches the 'a' values.
% pi_Smol_grid_ani_binary=zeros(size(Smol_grid_ani,1),size(Smol_grid_ani,1));
% % Lets also store the relevant z values while we are at it.
% pi_Smol_grid_ani_z_vals=struct();
% for ii=1:size(Smol_grid_ani,1) % This for-loop would be much simpler if I had stored the 'unique' indexes back when I created 'unique' from Smol_grid_ani_z
%     jj_c=1;
%     for jj=1:size(Smol_grid_ani,1)
%         % Smol_grid_ani has both 'a' and 'z'. If both of the 'a' are the
%         % same then we want to make this a 1.
%         if prod(Smol_grid_ani(ii,1:length(n_a))==Smol_grid_ani(jj,1:length(n_a)))==1 % If all of the 'a' elements are equal
%             pi_Smol_grid_ani_binary(ii,jj)=1;
%             if jj_c==1
%                 pi_Smol_grid_ani_z_vals(ii).val=Smol_grid_ani(jj,(length(n_a)+1):end);
%             else
%                 pi_Smol_grid_ani_z_vals(ii).val=[pi_Smol_grid_ani_z_vals(ii).val,Smol_grid_ani(jj,(length(n_a)+1):end)];                
%             end
%             jj_c=jj_c+1;
%         end
%     end
% end % This implementation is highly inefficient. No need for them to be nested in this way. (will only run once so can't be bothered fixing)
% % Now that we have all of these it is just a matter of going through each
% % of the "1's" and calculating the appropriate probability.
% for ii=1:size(Smol_grid_ani,1) % This for-loop would be much simpler if I had stored the 'unique' indexes back when I created 'unique' from Smol_grid_ani_z
%     z_ii=Smol_grid_ani(ii,(length(n_a)+1):end);
%     
%     clear z_grid_ani_ii % size changes in ii
%     z_grid_ani_ii=pi_Smol_grid_ani_z_vals(ii).val;
%     [z_grid_ani_ii_sorted,z_grid_ani_ii_sorted_ind]=sort(z_grid_ani_ii);
%     z_grid_ani_ii_sorted=(z_grid_ani_ii_sorted*(z_grid(2)-z_grid(1))+z_grid(1)+z_grid(2))/2; % Move off of hypercube.
%     
%     pi_z_ani_ii_sorted=nan(1,length(z_grid_ani_ii_sorted));
%     if length(pi_z_ani_ii_sorted)==1
%         pi_z_ani_ii_sorted=1;
%     else
%         % Use 'tauchen' method to calculate the transition probabilities
%         omega1=z_grid_ani_ii_sorted(2)-z_grid_ani_ii_sorted(1);
%         pi_z_ani_ii_sorted(1)=normcdf(z_grid_ani_ii_sorted(1)+omega1/2-rho*z_ii,mew,sigma);
%         for jj=2:length(z_grid_ani_ii_sorted)-1
%             omega1=z_grid_ani_ii_sorted(jj+1)-z_grid_ani_ii_sorted(jj);
%             omega2=z_grid_ani_ii_sorted(jj)-z_grid_ani_ii_sorted(jj-1);
%             pi_z_ani_ii_sorted(jj)=normcdf(z_grid_ani_ii_sorted(jj)+omega1/2-rho*z_ii,mew,sigma)-normcdf(z_grid_ani_ii_sorted(jj)-omega2/2-rho*z_ii,mew,sigma);
%         end
%         omega2=z_grid_ani_ii_sorted(end)-z_grid_ani_ii_sorted(end-1);
%         pi_z_ani_ii_sorted(end)=1-normcdf(z_grid_ani_ii_sorted(end)-omega2/2-rho*z_ii,mew,sigma);
%     end
%     
%     pi_z_ani_ii=pi_z_ani_ii_sorted(z_grid_ani_ii_sorted_ind);
%     jj_c=1;
%     for jj=1:size(Smol_grid_ani,1)
%         % Smol_grid_ani has both 'a' and 'z'. If both of the 'a' are the
%         % same then we want to make this a 1.
%         if pi_Smol_grid_ani_binary(ii,jj)==1 % If all of the 'a' elements are equal
%             pi_Smol_grid_ani(ii,jj)=pi_z_ani_ii(jj_c);
%             jj_c=jj_c+1;
%         end
%     end
%     
% end % This implementation is highly inefficient. No need for them to be nested in this way. (will only run once so can't be bothered fixing)


%% Now switch initial guess for value function (V_Dynare: based on 2nd order perturbation from Dynare) into the form of Smolyak-Chebyshev coefficients.

% Smol_grid_ani is currently in the [-1,1] hypercube. Move it to the full (a,z)-cube.
az_Smol_grid_ani = zeros(size(Smol_grid_ani));
for ii=1:length(n_a)
    amax=a_grid(sum(n_a(1:ii)));
    amin=a_grid(sum(n_a(1:ii))-1); % We know that n_a=2.
    az_Smol_grid_ani(:,ii)=(Smol_grid_ani(:,ii)*(amax-amin)+amin+amax)/2;
end
for ii=1:length(n_z)
    zmax=z_grid(sum(n_z(1:ii)));
    zmin=z_grid(sum(n_z(1:ii))-1); % We know that n_a=2.
    az_Smol_grid_ani(:,ii+length(n_a))=(Smol_grid_ani(:,ii+length(n_a))*(zmax-zmin)+zmin+zmax)/2;
end
% Numerical rounding error can leave the end points of az_Smol_grid_ani outside [amin,amax]x[zmin,zmax] by less than 10^(-16) 
% causing problems for interpolation later on. Following just corrects this when it occours.
for jj=1:size(az_Smol_grid_ani,1)
    for ii=1:length(n_a)
        amax=a_grid(sum(n_a(1:ii)));
        amin=a_grid(sum(n_a(1:ii))-1);
        if abs(az_Smol_grid_ani(jj,ii+length(n_a))-zmin)<10^(-15)
            az_Smol_grid_ani(jj,ii+length(n_a))=zmin;
        elseif abs(az_Smol_grid_ani(:,ii+length(n_a))-zmax)<10^(-15)
            az_Smol_grid_ani(jj,ii+length(n_a))=zmax;
        end
    end
    for ii=1:length(n_z)
        zmax=z_grid(sum(n_z(1:ii)));
        zmin=z_grid(sum(n_z(1:ii))-1);
        if abs(az_Smol_grid_ani(jj,ii+length(n_a))-zmin)<10^(-15)
            az_Smol_grid_ani(jj,ii+length(n_a))=zmin;
        elseif abs(az_Smol_grid_ani(:,ii+length(n_a))-zmax)<10^(-15)
            az_Smol_grid_ani(jj,ii+length(n_a))=zmax;
        end
    end
end


f_V0 = zeros(size(Smol_grid_ani,1),1);
[z_meshgrid, k_meshgrid]=meshgrid(z_grid_tauchen,k_gridtemp);
for ii=1:size(Smol_grid_ani,1)
    kz_sub=az_Smol_grid_ani(ii,:);
    f_V0(ii)=interp2(z_meshgrid, k_meshgrid,V_Dynare,kz_sub(2),kz_sub(1)); % Because meshgrid() works in a very stupid way in terms of dimensions (it thinks of column as the first dimension and row as the second; exactly the opposite of how they get indexed) this has to have kz_sub in the reversed ordering
end
% Evaluate the function on the Smolyak grid
Smol_polynom = Smolyak_Polynomial(Smol_grid_ani,numdimensions,mu_max,Smol_elem_ani); 
V0 = Smol_polynom\f_V0; % V is kept as the coefficients.

save ./SavedOutputs/StochNeoClassicalGrowth_SmolyakChebyV0.mat V0 f_V0
% load ./SavedOutputs/StochNeoClassicalGrowth_SmolyakChebyV0.mat V0 f_V0

%% Draw (another) graph of the initial value function


[z2,z1] = meshgrid((z_grid_tauchen-zmin)/(zmax-zmin),(k_gridtemp-kmin)/(kmax-kmin));
%[z1,z2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 
    
for j=1:size(z1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([z1(:,j),z2(:,j)],numdimensions,mu_max,Smol_elem_ani)*V0;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

[x1,x2] = meshgrid(linspace(z_grid(1),z_grid(2),size(z1,2)), linspace(k_grid(1),k_grid(2),size(z2,1)));

figure(7)
surf(x1,x2,V_fitted)
title('Value Function: Smolyak fit of V\_{}Dynare')    

figure(8)
surf(x1,x2,V_fitted-V_Dynare)
title('Value Function Check: V\_{}Dynare minus Smolyak fit of V\_{}Dynare')    


%% things to make my life easier while working on implementation of the codes
vfoptions.verbose=1;
Parameters=Params;
V0Kron=V0; % rather than V0Kron=V0, as V0 is not declare till a few lines later than this
VKron=V0Kron;
pi_z=pi_Smol_grid_ani;

%% Solve
%Do the value function iteration. Returns both the value function itself,
%and the optimal policy function.
d_grid=0; %no d variable
n_d=0; %no d variable

vfoptions.howards=10
vfoptions.maxhowards=0

tic;
%V0=nan; % If you do not wish to specify an intial guess just input V0=nan;
[V, Policy]=ValueFnIter_Case1(V0, n_d,n_k,n_z,d_grid,k_grid,z_grid, pi_z_unique, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
% pi_Smol_grid_ani
time=toc;

save ./SavedOutputs/StochNeoClassicalGrowth_SmolyakCheby.mat V Policy
% load ./SavedOutputs/StochNeoClassicalGrowth_SmolyakCheby.mat V Policy

fprintf('Time to solve the value function iteration was %8.2f seconds. \n', time)

    
%% Draw a graph of the value function
[z1,z2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 
    
for j=1:size(z1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([z1(:,j),z2(:,j)],numdimensions,mu_max,Smol_elem_ani)*V;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

[x1,x2] = meshgrid(linspace(z_grid(1),z_grid(2),size(z1,2)), linspace(k_grid(1),k_grid(2),size(z2,1)));

figure(2)
surf(x1,x2,V_fitted),title('Value Function: Smolyak interpolation')    
    
%% Draw a graph of the policy functions

for j=1:size(x1(1,:),2)
    Policy_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*Policy;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(4)
surf(x1,x2,Policy_fitted),title('Policy Function: Smolyak interpolation')

%% Draw graph of the value and policy functions as computed by VFI Toolkit
load ./GPU_Solution/Solution_gpu.mat V_gpu Policy_gpu k_grid_gpu z_grid_gpu

% Plot value function on the grid that is used by ChebySmol codes
figure(5)
surf(k_grid_gpu*ones(1,length(z_grid_gpu)),ones(length(k_grid_gpu),1)*z_grid_gpu',V_gpu)

% Plot policy function on the grid that is used by ChebySmol codes
figure(6)
surf(k_grid_gpu*ones(1,length(z_grid_gpu)),ones(length(k_grid_gpu),1)*z_grid_gpu',Policy_gpu)

    
%% Calculate Euler Equation Residuals




