%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'StochasticNeoclassicalGrowthModel';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('StochasticNeoclassicalGrowthModel.log');
M_.exo_names = 'epsilon';
M_.exo_names_tex = 'epsilon';
M_.exo_names_long = 'epsilon';
M_.endo_names = 'c';
M_.endo_names_tex = 'c';
M_.endo_names_long = 'c';
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.param_names = 'beta';
M_.param_names_tex = 'beta';
M_.param_names_long = 'beta';
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'gammac');
M_.param_names_tex = char(M_.param_names_tex, 'gammac');
M_.param_names_long = char(M_.param_names_long, 'gammac');
M_.param_names = char(M_.param_names, 'zbar');
M_.param_names_tex = char(M_.param_names_tex, 'zbar');
M_.param_names_long = char(M_.param_names_long, 'zbar');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'sigma_epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_epsilon');
M_.param_names_long = char(M_.param_names_long, 'sigma_epsilon');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 3;
M_.param_nbr = 7;
M_.orig_endo_nbr = 3;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('StochasticNeoclassicalGrowthModel_static');
erase_compiled_function('StochasticNeoclassicalGrowthModel_dynamic');
M_.lead_lag_incidence = [
 0 3 6;
 1 4 0;
 2 5 7;]';
M_.nstatic = 0;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 2;
M_.nspred   = 2;
M_.ndynamic   = 3;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(3, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(7, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 11;
M_.NNZDerivatives(2) = 14;
M_.NNZDerivatives(3) = -1;
load paramsfordynare;
set_param_value('beta',beta);
set_param_value('alpha',alpha);
set_param_value('delta',delta);
set_param_value('gammac',gammac);
set_param_value('rho',rho);
set_param_value('sigma_epsilon',sigma_epsilon);
M_.params( 5 ) = 0;
zbar = M_.params( 5 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(7)^2;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = 3.5;
oo_.steady_state( 2 ) = 50;
oo_.steady_state( 3 ) = 0;
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
options_.solve_algo = 2;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 50;
options_.order = 2;
options_.periods = 2000;
var_list_=[];
info = stoch_simul(var_list_);
save('StochasticNeoclassicalGrowthModel_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('StochasticNeoclassicalGrowthModel_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('StochasticNeoclassicalGrowthModel_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('StochasticNeoclassicalGrowthModel_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('StochasticNeoclassicalGrowthModel_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
