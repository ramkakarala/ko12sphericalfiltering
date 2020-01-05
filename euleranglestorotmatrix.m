function R = euleranglestorotmatrix(alpha,beta,gamma)
% usage
%        R = euleranglestorotmatrix(alpha,beta,gamma)
% r kakarala
%
% Provides the 3x3 rotation matrix R that provides the rotation
% defined by the Euler angles. Note: it is assumed that R is applied
% to row vectors x, i.e., xR = y.  The euler angles are 
% alpha = rotation around z axis
% beta = rotation around x axis following alpha
% gamma = rotation around z axis following alpha, beta

Ra = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];
Rb = [1 0 0; 0 cos(beta) sin(beta); 0 -sin(beta) cos(beta)];
Rg = [cos(gamma) sin(gamma) 0; -sin(gamma) cos(gamma) 0; 0 0 1];
R  = Ra*Rb*Rg;

