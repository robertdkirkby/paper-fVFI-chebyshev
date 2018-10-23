% Sobol (quasi-) random number sequence

% When doing high dimensional integration problems one option is
% Monte-Carlo. But rather than actually drawing from a random number
% generator we can use the quasi-random Sobol random number sequence which
% represents draws from d-dimensional Uniform(1,1) random distribution. By
% using Sobol we get the advantage that the quasi-random sequence has known 
% properties and so will typically
% acheive convergence of the integral faster than standard monte-carlo
% drawing from the same distribution would do.

% In two dimensions, first create a 'sobol set' object.
p=sobolset(2);
% Now get the first ten numbers from each (each is two-dimensional)
x=net(p,10);
