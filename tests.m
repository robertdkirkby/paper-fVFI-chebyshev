% For trying to figure out what is up with expections

size(pi_Smol_grid_ani) % Is the correct size
plot(sum(pi_Smol_grid_ani,2)) % All rows correctly sum to one

surf(cumsum(pi_Smol_grid_ani,2)) % Can see that most transitions are zero. Seem likely that this is the problem.


%% The transition matrix on the 'Chebyshev z-values' that I created based on Tauchen method.
pi_z_unique

Smol_grid_ani_z_unique_ab

%%

% The object stored as the Smolyak-Chebyshev polynomial coefficients is
VKron
% The object stored as values on the associated grid is
f_VKron

% My job is to calculate
E_f_VKron

% Currently my attempt is 
% E_f_VKron = pi_Smol_grid_ani*f_VKronOld; % Take expectation of values

% But this gives a 'wavy' version that is clearly not nearly accurate enough.

% 
z_grid_tauchen_hc=z_grid_tauchen/(((z_grid(2)-z_grid(1))+z_grid(1)+z_grid(2))/2);
pi_z_tauchen

a_grid_hc=sort(unique(Smol_grid_ani(:,1)));

z_grid_smolcheb_hc=sort(unique(Smol_grid_ani(:,2)));

pi_z_unique
surf(cumsum(pi_z_unique,2))

surf(f_VKron_fitted)

figure(6)
surf(E_f_VKron_fitted)


%% Draw a graph of the initial value function: V0

[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

V_fitted=nan(size(x1));    
for j=1:size(x1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*V0;
    % Evaluate Smolyak interpolating polynomial on the grid    
end


figure(1)
surf(x1,x2,V_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation

%% Draw a graph of the value function: VKron

[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

for j=1:size(x1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*VKron;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(2)
surf(x1,x2,V_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation

%% Draw a graph of the value function: VKronOld

[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

for j=1:size(x1(1,:),2)
    V_fittedOld(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*VKronOld;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(3)
surf(x1,x2,V_fittedOld),title('Smolyak interpolation')
    % Plot Smolyak interpolation

    
%% Draw a graph of the value function: f_VKron
VKronTemp = Smol_polynom\f_VKron; % VKron is outputted as the Smolyak-Chebyshev coefficients.

[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

for j=1:size(x1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*VKronTemp;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(2)
surf(x1,x2,V_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation

%% Draw a graph of the value function: f_VKron_fitted and E_f_VKron_fitted
[x1,x2] = meshgrid(a_grid_hc, z_grid_smolcheb_hc'); 
figure(6)
subplot(1,2,1); surf(x1,x2,f_VKron_fitted),title('f_VKron_fitted')
    % Plot Smolyak interpolation

subplot(1,2,2); surf(x1',x2',E_f_VKron_fitted),title('E_f_VKron_fitted')
    % Plot Smolyak interpolation


%% Draw a graph of the value function: E_f_VKron
VKronTemp3 = Smol_polynom\E_f_VKron; % VKron is outputted as the Smolyak-Chebyshev coefficients.

[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

for j=1:size(x1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*VKronTemp3;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(2)
surf(x1,x2,V_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation


%% Draw a graph of the value function: E_VKron

[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1] 

clear V_fitted
for j=1:size(x1(1,:),2)
    V_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],numdimensions,mu_max,Smol_elem_ani)*VKron3;
    % Evaluate Smolyak interpolating polynomial on the grid    
end

figure(3)
surf(x1,x2,V_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation

    
%% Presumably the error lies in pi_Smol_grid_ani

% Currently I do:
% E_f_VKron = pi_Smol_grid_ani*f_VKronOld; % Take expectation of values

% Lets start by just evaluating f_VKron and then doing a more brute force
% evaluation of E_fVKron from that.

























