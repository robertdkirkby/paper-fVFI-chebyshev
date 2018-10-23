% Neoclassical Stochastic Growth Model
% Example based on Diaz-Gimenez (2001) - Linear Quadratic Approximations: An Introduction 
% (Chapter 2 in 'Computational Methods for the Study of Dynamic Economies', edited by Marimon & Scott)

% This model is also used by Aldrich, Fernandez-Villaverde, Gallant, & Rubio-Ramirez (2011) - "Tapping the supercomputer under your desk: Solving dynamic equilibrium models with graphics processors,"
% But they do use slightly different parameters to those used here.

Javier=0;   %If you set this parameter to 0 then the parameters will all be set to those used by Aldrich, Fernandez-Villaverde, Gallant, & Rubio-Ramirez (2011)

%% Set up
tauchenoptions.parallel=2; % Use GPU (is anyway the default option)
vfoptions.parallel=2; % Use GPU (is anyway the default option)

% The sizes of the grids
n_z=25;
n_k=2^11;

load ../DynareInitialValueFns/paramsfordynare.mat beta gammac alpha delta rho sigma_epsilon
Params.beta=beta;
Params.gamma=gammac;
Params.alpha=alpha;
Params.delta=delta;
Params.rho=rho;
Params.sigma_epsilon=sigma_epsilon;

Params.sigmasq_epsilon=Params.sigma_epsilon^2;

%% Compute the steady state
K_ss=((alpha*beta)/(1-beta*(1-delta)))^(1/(1-alpha));
X_ss= delta*K_ss;
%These are not really needed; we just use them to determine the grid on
%capital. I mainly calculate them to stay true to original article.

%% Create grids (grids are defined as a column vectors)

q=4; % A parameter needed for the Tauchen Method
[z_grid, pi_z]=TauchenMethod(0,Params.sigmasq_epsilon,Params.rho,n_z,q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q), transmatix is (z,zprime)

k_grid=linspace(0,3*K_ss,n_k)'; % Grids should always be declared as column vectors

%% Now, create the return function
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val, s_val, gamma, alpha, delta) StochasticNeoClassicalGrowthModel_ReturnFn(aprime_val, a_val, s_val, gamma, alpha, delta);
ReturnFnParamNames={'gamma', 'alpha', 'delta'}; %It is important that these are in same order as they appear in 'StochasticNeoClassicalGrowthModel_ReturnFn'

%% Solve
%Do the value function iteration. Returns both the value function itself,
%and the optimal policy function.
d_grid=0; %no d variable
n_d=0; %no d variable

tic;
V0=ones(n_k,n_z);
[V, Policy]=ValueFnIter_Case1(V0, n_d,n_k,n_z,d_grid,k_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
time=toc;

fprintf('Time to solve the value function iteration was %8.2f seconds. \n', time)


%% Draw a graph of the value function

% Get [kmin,kmax] and [zmin,zmax] being used by ChebySmol codes
load ./gridendpoints.mat kmin kmax k_gridtemp zmin zmax z_grid_tauchen

[k_grid(1), k_grid(end); kmin kmax]
[z_grid(1), z_grid(end); zmin zmax]

% Graph of the value function on the grids that get used by the ChebySmol codes.
kminind=find(k_grid>kmin,1,'first')-1;
kmaxind=find(k_grid>kmax,1,'first');

k_grid_gpu=k_grid(kminind:kmaxind);
z_grid_gpu=z_grid;

% Plot value function on the grid that is used by ChebySmol codes
V_gpu=V(kminind:kmaxind,:);
surf(k_grid_gpu*ones(1,n_z),ones(length(k_grid_gpu),1)*z_grid_gpu',V_gpu)

% Plot policy function on the grid that is used by ChebySmol codes
Policy_gpu=shiftdim(Policy(1,kminind:kmaxind,:),1);
for ii=1:numel(Policy_gpu)
    Policy_gpu(ii)=k_grid(Policy_gpu(ii));
end
surf(k_grid_gpu*ones(1,n_z),ones(length(k_grid_gpu),1)*z_grid_gpu',Policy_gpu)

save ./Solution_gpu.mat V_gpu Policy_gpu k_grid_gpu z_grid_gpu

% % Or to get 'nicer' x and y axes use
% surf(k_grid*ones(1,n_z),ones(n_k,1)*z_grid',V)




