%% Using Smolyak-Chebyshev to fit the Value Function itself. What goes wrong

% Before running this you must run the following which creates V_Dynare and
% Policy_Dynare (and in turn requires that you run a Dynare mod-file)
% StochasticNeoClassicalGrowthModel_SmolyakChebyshev2.m

%% 1. Smolyak anisotropic method for 2 dimensions: set-up
% -----------------------------------------------
vector_mus_dimensions = [11,11]; % Introduce the level of approximation in every dimension from 1 to 10; see Section 4 of JMMV (2014)

d=length(vector_mus_dimensions); % Number of dimensions
mu_max  = max(vector_mus_dimensions); % Compute the maximum level of approximation across all dimensions

Smolyak_elem_iso = Smolyak_Elem_Isotrop(d,mu_max);
    % Construct the matrix of indices of multidimesional Smolyak elements (grid points and polynomial basis functions) that satisfy the usual
    % isotropic Smolyak rule for the approximation level equal to "mu_max"
Smol_elem_ani = Smolyak_Elem_Anisotrop(Smolyak_elem_iso,vector_mus_dimensions);
    % Select from the matrix of isotropic indices "Smol elem_iso" a subset of indices that correspond to the given anisotropic "vector_mus_dimensions"
Smol_grid_ani = Smolyak_Grid(d,mu_max,Smol_elem_ani);
    % Construct the Smolyak grid for the given subindices of anisotropic Smolyak elements "Smol_elem_ani"
Smol_polynom = Smolyak_Polynomial(Smol_grid_ani,d,mu_max,Smol_elem_ani);
    % Matrix of the polynomial basis functions evaluated in the grid points

if d==2
    figure(6)
    scatter(Smol_grid_ani(:,1),Smol_grid_ani(:,2),'filled'),title('Smolyak grid') % Enable it to draw the Smolyak grid
end

%% 2. Complicated Example: Interpolation of a function y=2*x1 .*exp(-4*x1.^2-16*x2.^2);
% ---------------------------------------------------------------
y_Smolyak = 2*Smol_grid_ani(:,1) .* exp(-4*Smol_grid_ani(:,1).^2 - 16*Smol_grid_ani(:,2).^2);
    % Evaluate the function on the Smolyak grid

b = Smol_polynom\y_Smolyak;
    % Compute the coefficients of Smolyak interpolating polynomial

size(b)
    % Take a look to see how many coefficients this involves

% Compare the true and interpolated functions on a dense grid
% ---------------------------------------------------------------
[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1]

y_true = 2*x1.* exp(-4*x1.^2 - 16*x2.^2);
    % Evaluate the true function on the grid

figure(1)
subplot(1,2,1), surf(x1,x2,y_true),title('True function')
    % Plot the true function

clear y_fitted
for j=1:size(x1(1,:),2)
    y_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],d,mu_max,Smol_elem_ani)*b;
    % Evaluate Smolyak interpolating polynomial on the grid
end

subplot(1,2,2), surf(x1,x2,y_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation

LInf=max(max(abs(y_true-y_fitted)))

% figure(12)
% subplot(2,2,1), surf(x1,x2,y_true, 'EdgeColor','none'),title('True function')
% subplot(2,2,2), surf(x1,x2,y_fitted22, 'EdgeColor','none'),title('Smolyak-Chebyshev Approximation: 2,2')
% subplot(2,2,3), surf(x1,x2,y_fitted55, 'EdgeColor','none'),title('Smolyak-Chebyshev Approximation: 5,5')
% subplot(2,2,4), surf(x1,x2,y_fitted77, 'EdgeColor','none'),title('Smolyak-Chebyshev Approximation: 7,7')


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

figure(2)
subplot(1,2,1), surf(x1,x2,y_true),title('True function')
    % Plot the true function

clear y_fitted
for j=1:size(x1(1,:),2)
    y_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],d,mu_max,Smol_elem_ani)*b;
    % Evaluate Smolyak interpolating polynomial on the grid
end

subplot(1,2,2), surf(x1,x2,y_fitted),title('Smolyak interpolation')

LInf=max(max(abs(y_true-y_fitted)))

y_error=y_true-y_fitted;

figure(5)
surf(x1,x2,y_error)
title('True minus Smolyak')
    % Plot Smolyak interpolation
    
% Figures for paper
figure(10)
surf(x1,x2,y_true)
figure(11)
surf(x1,x2,y_fitted)

% figure(12)
% subplot(2,2,1), surf(x1,x2,y_true, 'EdgeColor','none'),title('Value function')
% subplot(2,2,2), surf(x1,x2,y_error22, 'EdgeColor','none'),title('Approximation Error: Order 2,2')
% subplot(2,2,3), surf(x1,x2,y_error55, 'EdgeColor','none'),title('Approximation Error: Order 5,5')
% subplot(2,2,4), surf(x1,x2,y_error77, 'EdgeColor','none'),title('Approximation Error: Order 7,7')

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

figure(3)
subplot(1,2,1), surf(x1,x2,y_true),title('True function')
    % Plot the true function

clear y_fitted
for j=1:size(x1(1,:),2)
    y_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],d,mu_max,Smol_elem_ani)*b;
    % Evaluate Smolyak interpolating polynomial on the grid
end

subplot(1,2,2), surf(x1,x2,y_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation
    
LInf=max(max(abs(y_true-y_fitted)))
    
figure(6)
y_error=y_true-y_fitted;
surf(x1,x2,y_error)
title('True minus Smolyak')


% figure(12)
% subplot(2,2,1), surf(x1,x2,y_true, 'EdgeColor','none'),title('Policy function')
% subplot(2,2,2), surf(x1,x2,y_error22, 'EdgeColor','none'),title('Approximation Error: Order 2,2')
% subplot(2,2,3), surf(x1,x2,y_error55, 'EdgeColor','none'),title('Approximation Error: Order 5,5')
% subplot(2,2,4), surf(x1,x2,y_error77, 'EdgeColor','none'),title('Approximation Error: Order 7,7')
