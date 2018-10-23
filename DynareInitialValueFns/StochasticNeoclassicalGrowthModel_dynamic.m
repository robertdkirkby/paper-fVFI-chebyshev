function [residual, g1, g2, g3] = StochasticNeoclassicalGrowthModel_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(3, 1);
T15 = params(1)*y(6)^(-params(4));
T23 = params(3)*exp(y(7))*y(4)^(params(3)-1);
T26 = 1+T23-params(2);
T47 = params(1)*getPowerDeriv(y(6),(-params(4)),1);
lhs =y(3)^(-params(4));
rhs =T15*T26;
residual(1)= lhs-rhs;
lhs =y(3)+y(4);
rhs =exp(y(5))*y(1)^params(3)+y(1)*(1-params(2));
residual(2)= lhs-rhs;
lhs =y(5);
rhs =params(6)*y(2)+x(it_, 1);
residual(3)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(3, 8);

  %
  % Jacobian matrix
  %

  g1(1,3)=getPowerDeriv(y(3),(-params(4)),1);
  g1(1,6)=(-(T26*T47));
  g1(1,4)=(-(T15*params(3)*exp(y(7))*getPowerDeriv(y(4),params(3)-1,1)));
  g1(1,7)=(-(T15*T23));
  g1(2,3)=1;
  g1(2,1)=(-(1-params(2)+exp(y(5))*getPowerDeriv(y(1),params(3),1)));
  g1(2,4)=1;
  g1(2,5)=(-(exp(y(5))*y(1)^params(3)));
  g1(3,2)=(-params(6));
  g1(3,5)=1;
  g1(3,8)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  v2 = zeros(14,3);
  v2(1,1)=1;
  v2(1,2)=19;
  v2(1,3)=getPowerDeriv(y(3),(-params(4)),2);
  v2(2,1)=1;
  v2(2,2)=46;
  v2(2,3)=(-(T26*params(1)*getPowerDeriv(y(6),(-params(4)),2)));
  v2(3,1)=1;
  v2(3,2)=30;
  v2(3,3)=(-(T47*params(3)*exp(y(7))*getPowerDeriv(y(4),params(3)-1,1)));
  v2(4,1)=1;
  v2(4,2)=44;
  v2(4,3)=  v2(3,3);
  v2(5,1)=1;
  v2(5,2)=28;
  v2(5,3)=(-(T15*params(3)*exp(y(7))*getPowerDeriv(y(4),params(3)-1,2)));
  v2(6,1)=1;
  v2(6,2)=54;
  v2(6,3)=(-(T23*T47));
  v2(7,1)=1;
  v2(7,2)=47;
  v2(7,3)=  v2(6,3);
  v2(8,1)=1;
  v2(8,2)=52;
  v2(8,3)=(-(T15*params(3)*exp(y(7))*getPowerDeriv(y(4),params(3)-1,1)));
  v2(9,1)=1;
  v2(9,2)=31;
  v2(9,3)=  v2(8,3);
  v2(10,1)=1;
  v2(10,2)=55;
  v2(10,3)=(-(T15*T23));
  v2(11,1)=2;
  v2(11,2)=1;
  v2(11,3)=(-(exp(y(5))*getPowerDeriv(y(1),params(3),2)));
  v2(12,1)=2;
  v2(12,2)=33;
  v2(12,3)=(-(exp(y(5))*getPowerDeriv(y(1),params(3),1)));
  v2(13,1)=2;
  v2(13,2)=5;
  v2(13,3)=  v2(12,3);
  v2(14,1)=2;
  v2(14,2)=37;
  v2(14,3)=(-(exp(y(5))*y(1)^params(3)));
  g2 = sparse(v2(:,1),v2(:,2),v2(:,3),3,64);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],3,512);
end
end
