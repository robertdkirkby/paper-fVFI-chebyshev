function [L2_V,LInf_V, L2_Policy, LInf_Policy]=WhatGoesWrong_Table_fn(vector_mus_dimensions, use_smolyak)
% Pretty much same as 'WhatGoesWrong.m', but now takes
% vector_mus_dimensions as input, and just returns some measures used for a
% table.

%% 1. Smolyak anisotropic method for 2 dimensions: set-up
% -----------------------------------------------

% This is now provided as function input
% vector_mus_dimensions = [2,2]; % Introduce the level of approximation in every dimension from 1 to 10; see Section 4 of JMMV (2014)

d=length(vector_mus_dimensions); % Number of dimensions
if use_smolyak==0
    mu_max  = max(vector_mus_dimensions); % Compute the maximum level of approximation across all dimensions
elseif use_smolyak==1
    mu_max = ceil(sqrt(max(max(vector_mus_dimensions))));
end

Smolyak_elem_iso = Smolyak_Elem_Isotrop(d,mu_max);
    % Construct the matrix of indices of multidimesional Smolyak elements (grid points and polynomial basis functions) that satisfy the usual
    % isotropic Smolyak rule for the approximation level equal to "mu_max"
Smol_elem_ani = Smolyak_Elem_Anisotrop(Smolyak_elem_iso,vector_mus_dimensions);
    % Select from the matrix of isotropic indices "Smol elem_iso" a subset of indices that correspond to the given anisotropic "vector_mus_dimensions"
Smol_grid_ani = Smolyak_Grid(d,mu_max,Smol_elem_ani);
    % Construct the Smolyak grid for the given subindices of anisotropic Smolyak elements "Smol_elem_ani"
Smol_polynom = Smolyak_Polynomial(Smol_grid_ani,d,mu_max,Smol_elem_ani);
    % Matrix of the polynomial basis functions evaluated in the grid points
    
%% 3. Example: Fit Value Function for Stochastic Neoclassical Growth Model
% ---------------------------------------------------------------
load ./SavedOutput/StochNeoClassicalGrowth_VfromDynare.mat V_Dynare k_gridtemp z_grid_tauchen

y_V_Dynare = zeros(size(Smol_grid_ani,1),1);
zmin=min(z_grid_tauchen); zmax=max(z_grid_tauchen); kmin=min(k_gridtemp); kmax=max(k_gridtemp);
[z_meshgrid_hc, k_meshgrid_hc]=meshgrid((2*z_grid_tauchen-zmax-zmin)/(zmax-zmin),(2*k_gridtemp-kmax-kmin)/(kmax-kmin));
for ii=1:size(Smol_grid_ani,1)
    k_val_hc=Smol_grid_ani(ii,1);
    z_val_hc=Smol_grid_ani(ii,2);
    y_V_Dynare(ii)=interp2(z_meshgrid_hc, k_meshgrid_hc,V_Dynare,z_val_hc,k_val_hc); % Because meshgrid() works in a very stupid way in terms of dimensions (it thinks of column as the first dimension and row as the second; exactly the opposite of how they get indexed) this has to have kz_sub in the reversed ordering
end
    % Evaluate the function on the Smolyak grid

b = Smol_polynom\y_V_Dynare;
    % Compute the coefficients of Smolyak interpolating polynomial

size(b)
    % Take a look to see how many coefficients this involves

% Compare the true and interpolated functions on a dense grid
% ---------------------------------------------------------------
% [x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
x1=k_meshgrid_hc; x2=z_meshgrid_hc;
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1]

y_true = V_Dynare;
    % Evaluate the true function on the grid

clear y_fitted
for j=1:size(x1(1,:),2)
    y_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],d,mu_max,Smol_elem_ani)*b;
    % Evaluate Smolyak interpolating polynomial on the grid
end

y_approximationerrors=y_true-y_fitted;

L2_V=sqrt(mean(mean(y_approximationerrors.^2)));
LInf_V=max(max(abs(y_approximationerrors)));

%% 3. Example: Fit Policy Function for Stochastic Neoclassical Growth Model
% ---------------------------------------------------------------
load ./SavedOutput/StochNeoClassicalGrowth_VfromDynare.mat Policy_Dynare k_gridtemp z_grid_tauchen

y_Policy_Dynare = zeros(size(Smol_grid_ani,1),1);
zmin=min(z_grid_tauchen); zmax=max(z_grid_tauchen); kmin=min(k_gridtemp); kmax=max(k_gridtemp);
[z_meshgrid_hc, k_meshgrid_hc]=meshgrid((2*z_grid_tauchen-zmax-zmin)/(zmax-zmin),(2*k_gridtemp-kmax-kmin)/(kmax-kmin));
for ii=1:size(Smol_grid_ani,1)
    k_val_hc=Smol_grid_ani(ii,1);
    z_val_hc=Smol_grid_ani(ii,2);
    y_Policy_Dynare(ii)=interp2(z_meshgrid_hc, k_meshgrid_hc,Policy_Dynare,z_val_hc,k_val_hc); % Because meshgrid() works in a very stupid way in terms of dimensions (it thinks of column as the first dimension and row as the second; exactly the opposite of how they get indexed) this has to have kz_sub in the reversed ordering
end
    % Evaluate the function on the Smolyak grid

b = Smol_polynom\y_Policy_Dynare;
    % Compute the coefficients of Smolyak interpolating polynomial

size(b)
    % Take a look to see how many coefficients this involves

% Compare the true and interpolated functions on a dense grid
% ---------------------------------------------------------------
% [x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
x1=k_meshgrid_hc; x2=z_meshgrid_hc;
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1]

y_true = Policy_Dynare;
    % Evaluate the true function on the grid

clear y_fitted
for j=1:size(x1(1,:),2)
    y_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],d,mu_max,Smol_elem_ani)*b;
    % Evaluate Smolyak interpolating polynomial on the grid
end

y_approximationerrors=y_true-y_fitted;

L2_Policy=sqrt(mean(mean(y_approximationerrors.^2)));
LInf_Policy=max(max(abs(y_approximationerrors)));

