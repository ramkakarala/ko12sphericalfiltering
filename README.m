%% MATLAB software accompanying the paper
% "A phase-sensitive approach to filtering on the sphere"
% R Kakarala and P. Ogunbona
% Submitted to IEEE Trans Signal Processing
% For review purposes only
% Included in this directory are software written by 
% Professor Moo K. Chung, 
% http://www.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/weighted-SPHARM.html
% and
% Professor F. J. Simons
% http://geoweb.princeton.edu/people/simons/software.html
% Please refer to the above webpages for more details.
%% 
% The main script is "experimentstranssp.m" file.  Before running it, 
% be sure to run the "SPHARMconstruct.m" to create and store the spherical
% harmonic basis as follows
%%
SPHARMconstruct(65);
%% main file
experimentstranssp;