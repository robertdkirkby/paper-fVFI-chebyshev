%%% Example using a variant of the Basic RBC model (following Aruoba, Fernandez-Villaverde, & Rubio-Ramirez, 2006)

tic;
UseAlternativeParams=1;

%% Set up
tauchenoptions.parallel=0;
vfoptions.parallel=0;

q=3; %Parameter for the Tauchen method

%Discounting rate
Params.beta = 0.9896;

%Parameter values
Params.alpha = 0.4; % alpha
Params.theta=0.357;
Params.rho = 0.95; % rho
Params.tau=2;
Params.delta = 0.0196; % delta
Params.sigmasq_epsilon=0.007^2;

if UseAlternativeParams==1
    Params.tau=50;
    Params.sigmasq_epsilon=0.035^2;
end

%% Compute the steady state (just use these when picking grids)
varphi=(1/Params.alpha*(1/Params.beta-1+Params.delta))^(1/(1-Params.alpha));
Psi=Params.theta/(1-Params.theta)*(1-Params.alpha)*varphi^(-Params.alpha);
Omega=varphi^(1-Params.alpha)-Params.delta;
K_ss=Psi/(Omega+varphi*Psi);


%% Create grids (it is very important that each grid is defined as a column vector)

% Using Smolyak-Chebyshev method to approximate the value function and policy function.
vfoptions.solnmethod='smolyak_chebyshev';
% Set 'level of approximation' for each dimension of the value function (must be given for (d,a,z)).
vfoptions.aniso_levelofapprox=[7,3]; % If you don't set this it is assumed you want '5' for each and every dimension.

a_grid=[0.01; 5*K_ss];
d_grid=[0;1];

n_d=2;
n_a=2;

n_z=21;  %Note, just need one for tauchen method
[z_grid, ~]=TauchenMethod(0,sigmasq_epsilon,rho,n_z,q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
% Note that we won't actually use this Tauchen AR(1). But we will use it to get max and min z values.

z_grid=[z_grid(1),z_grid(end)]'; n_z=2;
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
    pi_z_unique(ii,end)=1-normcdf(Smol_grid_ani_z_unique_ab(jj)-omega2/2-rho*Smol_grid_ani_z_unique_ab(ii),mew,sigma);
end
% (given that Tauchen method is normally implemented with matrices not for loops this could certainly be improved)
% For more general case, maybe a simulation combined with nearest
% neighbour to assign each simulated observation with a point, then divide by number 
% by number of simulations (do for each starting point). This is not going
% to be the fastest way, but will work.
%
% Now that we have the pi_z_unique, we just create the full pi_Smol_grid_ani which
% gives the transition probabilities between the various rows of Smol_grid_ani and
% populate it using the pi_z_unique values.
pi_Smol_grid_ani=zeros(size(Smol_grid_ani,1),size(Smol_grid_ani,1));
for ii=1:size(Smol_grid_ani,1) % This for-loop would be much simpler if I had stored the 'unique' indexes back when I created 'unique' from Smol_grid_ani_z
    for ii2=1:size(Smol_grid_ani_z_unique,1)
        if Smol_grid_ani_z(ii,:)==Smol_grid_ani_z_unique(ii2,:)
            ii_c=ii2;
        end
    end
    for jj=1:size(Smol_grid_ani,1)
        for jj2=1:size(Smol_grid_ani_z_unique,1)
            if Smol_grid_ani_z(jj,:)==Smol_grid_ani_z_unique(jj2,:)
                jj_c=jj2;
            end
        end
        % Smol_grid_ani has both 'a' and 'z'. We have found the 'z' and
        % 'zprime' and we have the probability of moving between stored in
        % pi_z_unique. If both of the 'a' are the same then we want to
        % store this value in pi_Smol_grid_ani.
        if prod(Smol_grid_ani(ii,1:length(n_a))==Smol_grid_ani(jj,1:length(n_a)))==1 % If all of the 'a' elements are equal
            pi_Smol_grid_ani(ii,jj)=pi_z_unique(ii_c,jj_c);
        end
    end
end % This implementation is highly inefficient. No need for them to be nested in this way. (will only run once so can't be bothered fixing)
for ii=1:size(Smol_grid_ani,1)
    pi_Smol_grid_ani(ii,:)=pi_Smol_grid_ani(ii,:)/sum(pi_Smol_grid_ani(ii,:)); % Normalize rows so that each row sums to one.
end % This loop is unnecessary, could be done as matrix operation
%
% Phew! Done. Now have pi_z for a sparse Smolyak-Chebyshev grid!!!
% surf(cumsum(pi_Smol_grid_ani,2))


%% Now, create the return function 
DiscountFactorParamNames={'beta'};

ReturnFn=@(d_val,aprime_val, a_val, s_val,alpha,delta,theta,tau) BasicRealBusinessCycleModel_ReturnFn(d_val,aprime_val, a_val, s_val,alpha,delta,theta,tau);
ReturnFnParams={'alpha','delta','theta','tau'}; %It is important that these are in same order as they appear in 'BasicRealBusinessCycleModel_ReturnFn'

%% Things to make my life easier while working on implementation of the codes
vfoptions.verbose=1;
Parameters=Params;
V0Kron=nan; % rather than V0Kron=V0, as V0 is not declare till a few lines later than this
VKron=V0Kron;

%% Use Dynare to solve model and get an initial guess of value function.
% Dynare solution using second-order perturbation methods.

% Save the parameters so that dynare can use them to solve 2nd perturbation version of model which is used to get initial guess of value function
beta=Params.beta;
gammac=Params.gamma;
alpha=Params.alpha;
delta=Params.delta;
Params.theta;
Params.tau;
rho=Params.rho;
sigma_epsilon=Params.sigma_epsilon;
save ./DynareInitialValueFns/paramsfordynareRBC.mat beta gammac alpha delta rho sigma_epsilon

% Since the maintainers of Dynare are a bunch of useless bums who clearly
% don't use matlab and ubuntu the 'dynare-matlab' package is out of date
% and can't be installed. Instead you will just have to manually run dynare
% from Octave after having created 'paramsfordynare.mat', and then come
% back here to load up oo_. and get the policy functions from that.
% dynare ./DynareInitialValueFns/BasicRBC_SocialPlanner.mod

%% Solve
%Do the value function iteration. Returns both the value function itself, and the optimal policy function.
V0=nan;
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_k,n_z,d_grid,k_grid,z_grid, pi_Smol_grid_ani, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
time1=toc;

save ./SavedOutputs/BasicRBC_SmolyakCheby.mat V Policy

fprintf('Run time for value function iteration: %8.2f seconds \n', time1)


%% Draw a graph of the value function

[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

for j=1:size(x1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*V;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(1)
surf(x1,x2,V_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation


%% Draw a graph of the policy functions

for j=1:size(x1(1,:),2)
    aPolicy_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*Policy(:,1);
    dPolicy_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*Policy(:,2);
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(2)
sub(2,1,1); surf(x1,x2,aPolicy_fitted),title('aprime: Smolyak interpolation')
sub(2,1,2); surf(x1,x2,dPolicy_fitted),title('d: Smolyak interpolation')
    % Plot Smolyak interpolation

%% Calculate Euler Equation Residuals


















