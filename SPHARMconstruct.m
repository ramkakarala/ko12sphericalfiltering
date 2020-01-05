function []=SPHARMconstruct(L)
% usage
%  SPHARMconstruct(L)
%
%
%                 WARNINING: The hard drive must have at least 2.4GB of space.
%  L           : degree of spherical harmonics. The recommend value of k is less than 86 otherwise 
%                 you will encounter the numerical singularity caused by MATLAB. 
%
%  Consruct spherical harmonic functions upto the degree k. 
%  Requires two external files Y_l.m and unitsphere.mat
%
% (C) Moo K. Chung, 2006-2008
%       Department of Biostatistics and Medical Informatics
%       University of Wisconsin-Maison
%  
% email://mkchung@wisc.edu
% http://www.stat.wisc.edu/softwares/weighted-SPHARM/weighted-SPHARM.html
%
% If you use this code, please reference the following paper. 
% You need to read the paper to modify the code.
%
% Chung, M.K., Dalton, K.M., Shen, L., L., Evans, A.C., Davidson, R.J. 2007. 
% Weighted Fourier series representation and its application to quantifying 
% the amount of gray matter. IEEE Transactions on Medical Imaging, 26:566-581.
%
% Update history: 
% Created Sept 19 2006; Modified July 5, 2007 
%-------------------------------------------------------------------------------

load unitsphere.mat
% if directory(end)~='\'
%     directory(end+1)='\';
% end;
for l=0:L
    %elapsedtime(l+1)=toc; these lines for measuring running time
    %[l elapsedtime(l+1)]
    temp=Y_l(l,theta,varphi);
    
    f_name=int2str(l); %%strcat(directory,int2str(l));
    save(f_name,'temp')
    fprintf(1,'index l=%d done\n',l);
end;
fprintf(1,'All done\n');


