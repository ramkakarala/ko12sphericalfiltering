%% This file contains the matlab code used for the experiments described
% in the paper "A phase-sensitive approach to filtering on the sphere" 
% by R. Kakarala and P. Ogunbona, submitted to IEEE Transactions on 
% Signal Processing
% For review purposes only
% ramakrishna@ntu.edu.sg
%%
clear;
close all;
set(0,'defaulttextinterpreter','latex')
%% now read in the world map
F = readS2Kitfdata('mapofworld.dat');
BW = max(size(F))/2;
%% compute angular samples
[x,y,z]=sphere(2*BW-1);
Npts = 2*BW;
jkvec = 0:2*BW-1;
thetavec = pi/(4*BW)*(2*jkvec+1);
varphivec = 2*pi*jkvec/(2*BW);
[varphi,theta]=meshgrid([varphivec 2*pi],[0 thetavec pi]);
% note the addition of 0, pi, and 2pi above makes the spherical
% triangulation complete. Those points are NOT sampled by S2kit
% for purposes of computing Spherical harmonic transform
tri=delaunay(theta,varphi);
x=cos(varphi).*sin(theta);
y=sin(varphi).*sin(theta);
z=cos(theta);
%% extrapolate to fill in the values at 0 and pi radians (latitude)
Fe = zeros(size(theta));
Fe(2:2*BW+1,1:2*BW) = F;
% the last column of Fe is at varphi = 2pi, same as 0
Fe(:,end) = Fe(:,1);
% first row of Fe is where theta = 0; extrapolate it by adding deriv
linearexp = F(1,:)+F(1,:)-F(2,:);
Fe(1,:) = mean(linearexp)*ones(1,size(theta,2));
% last row of Fe is where theta = pi; extrapolate it
linearexp = F(end,:)+F(end,:)-F(end-1,:);
Fe(end,:) = mean(linearexp)*ones(1,size(theta,2));
%%
mF = min(Fe(:));
MF = max(Fe(:));
Fs = Fe - min(Fe(:))+1.0;  % now a function value of zero means it lies on the sphere
%Fs = 1 + (Fe-mF)/(MF-mF); % this is an alternative normalization
figure
%subplot(1,2,1);
trisurf(tri,x.*Fs,y.*Fs,z.*Fs,Fs);  % plot function on sphere
%title('World');
axis equal
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%% Fourier expansion of world map
%subplot(1,2,2);
L=BW-1;
clear M;
ll = L;

sigma=0.0005; %important value, chosen after some experiments 
% value determines the degree to which Gibbs phenomenon is suppressed

[Fsmoothworld,fourier_coeffsworld]= SPHARMsmoothscalarcomplex(Fe,ll,sigma,theta,varphi,'normal');
if (norm(imag(Fsmoothworld))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;

% plot result of fitting with L spherical harmonics

figure;
Fsmoothr=real(Fsmoothworld);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('World, %d-th order approximation',ll));
%% Fourier coeffs of impulse response
fourier_coeffsimpulse = zeros(size(fourier_coeffsworld));
L = BW-1;
for k=0:L-1
    fourier_coeffsimpulse(k+1,1)=sqrt((2*k+1)/(4*pi));
end;
fourier_coeffsimpulse_norm = fourier_coeffsimpulse/sum(fourier_coeffsimpulse(:));
fimpulserecon = SPHARMreconscalarcomplex(fourier_coeffsimpulse_norm,L,sigma,theta,varphi,'normal');
if (norm(imag(fimpulserecon))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
%subplot(1,2,1);
Fsmoothr=real(fimpulserecon);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF)/(MF-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%colormap bone
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('impulse function',ll));
%% create wigner matrices for rotation
L = BW-1; %input('maximum order L=');
beta = 0; alpha=pi/8; gamma=0; % rotation around the north pole
% check if the reps have the right character
DC = dmatrixbeta(L,beta);
charok = 1;
for k=0:L-1
    ind = -k:k;
    D = exp(j*ind'*alpha)*exp(j*ind*gamma) .* DC{k+1};
    t = real(trace(D));
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)
        ch = sin((k+0.5)*phi)/sin(phi/2);
    else
        if k==1
            phi = acos((t-1)/2);
            if (phi == 0)
                ch = 2*k+1;
            else
                ch = sin((k+0.5)*phi)/sin(phi/2);
            end;
        else
            ch = 1;
        end;
    end;
    if (abs (t-ch) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
    Dmat{k+1}=D;
end;
if (charok)
    fprintf(1,'Character test passed\n');
end;
%% rotate using Wigner matrices
Frot = zeros(size(fourier_coeffsworld));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsworld(k+1,1:flen);
    % coeffs are [0,+1,+2,..+ell,-1,-2...-ell]
    % must unwrap them into [-ell,..,-1,0,1,...ell];
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    frotvec = funwrap*Dmat{k+1};
    % now rewrap them to store;
    frotwrap = [frotvec(k+1),frotvec(k+2:end),fliplr(frotvec(1:k))];
    Frot(k+1,1:flen) = frotwrap;
end;
%% reconstruct the rotated function from coefficients only
frotrecon = SPHARMreconscalarcomplex(Frot,L,sigma,theta,varphi,'normal');
if (norm(imag(frotrecon))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
%subplot(2,2,3);
Fsmoothr=real(frotrecon);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
title(sprintf('rotated original, recon L=%d',ll));
%% create new wigner matrices for rotation
L = BW-1; %input('maximum order L=');
beta = pi/8; alpha=pi/2; gamma=-pi/2; % rotation around the east-west axis
% check if the reps have the right character
DC = dmatrixbeta(L,beta);
charok = 1;
for k=0:L-1
    ind = -k:k;
    D2 = exp(j*ind'*alpha)*exp(j*ind*gamma) .* DC{k+1};
    t = real(trace(D2));
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)
        ch = sin((k+0.5)*phi)/sin(phi/2);
    else
        if k==1
            phi = acos((t-1)/2);
            if (phi == 0)
                ch = 2*k+1;
            else
                ch = sin((k+0.5)*phi)/sin(phi/2);
            end;
        else
            ch = 1;
        end;
    end;
    if (abs (t-ch) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
    Dmat2{k+1}=D2;
end;
if (charok)
    fprintf(1,'Character test passed\n');
end;
%% rotate using new Wigner matrices
Frot2 = zeros(size(fourier_coeffsworld));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsworld(k+1,1:flen);
    % coeffs are [0,+1,+2,..+ell,-1,-2...-ell]
    % must unwrap them into [-ell,..,-1,0,1,...ell];
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    frotvec = funwrap*Dmat2{k+1};
    % now rewrap them to store;
    frotwrap = [frotvec(k+1),frotvec(k+2:end),fliplr(frotvec(1:k))];
    Frot2(k+1,1:flen) = frotwrap;
end;
%% reconstruct the new rotated function from coefficients only
frotrecon2 = SPHARMreconscalarcomplex(Frot2,L,sigma,theta,varphi,'normal');
if (norm(imag(frotrecon2))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
%subplot(2,2,4);
Fsmoothr=real(frotrecon2);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
title(sprintf('rotated original with new rotation, recon L=%d',ll));
%% create simple 5 pt filter both lat and long lines filter
% alpha is rotation around z axis, i.e., along lines of latitude
% beta is rotation around new x axis
% gamma is rotation around new z axis
beta1 = pi/32; alpha1=pi/2; gamma1=-pi/2;  % latitude only
beta2 = pi/32; alpha2 = 0; gamma2=0; % longitude only
%% check if the reps have the right character
DC1 = dmatrixbeta(L,beta1);  % identify
DC2 = dmatrixbeta(L,beta2);
charok = 1;
for k=0:L-1
    ind = -k:k;
    D1 = exp(j*ind'*alpha1)*exp(j*ind*gamma1) .* DC1{k+1};
    D2 = exp(j*ind'*alpha2)*exp(j*ind*gamma2) .* DC2{k+1};
    
    t1 = real(trace(D1));
    
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)
        ch = sin((k+0.5)*phi)/sin(phi/2);
    else
        if k==1
            phi = acos((t1-1)/2);
            if (phi == 0)
                ch = 2*k+1;
            else
                ch = sin((k+0.5)*phi)/sin(phi/2);
            end;
        else
            ch = 1;
        end;
    end;
    if (abs (t1-ch) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
    t2 = real(trace(D2));
    
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)
        ch2 = sin((k+0.5)*phi2)/sin(phi2/2);
    else
        if k==1
            phi2 = acos((t2-1)/2);
            if (phi2 == 0)
                ch2 = 2*k+1;
            else
                ch2 = sin((k+0.5)*phi2)/sin(phi2/2);
            end;
        else
            ch2 = 1;
        end;
    end;
    if (abs (t2-ch2) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
    Dfiltmat4pt{k+1}=(4*eye(size(D1))+D1+D1'+D2+D2')/8;
end;
if (charok)
    fprintf(1,'Character test passed\n');
end;
%% filter's impulse response
Ffilt4pt = zeros(size(fourier_coeffsimpulse));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsimpulse(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Ffilt4pt(k+1,1:flen) = ffiltwrap;
end;
%% recon filtered
ffiltrecon4pt = SPHARMreconscalarcomplex(Ffilt4pt,L,sigma,theta,varphi,'normal');
if (norm(imag(ffiltrecon4pt))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
Fsmoothr=real(ffiltrecon4pt);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF)/(MF-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
%camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%% plot magnitude response of laplacian
maglap4pt=zeros(1,L);
for k=0:L-1
    maglap4pt(k+1) = norm(Ffilt4pt(k+1,:))/sqrt((2*k+1)/(4*pi));
    % division normalizes by magnitude of impulse function
end;
figure,plot(maglap4pt),xlabel('$\ell$'),ylabel('$\| F \|$');
%title(sprintf('filtered along 4pt L=%d',ll));
%% Filter world map with 4pt
Fworldfilt4pt = zeros(size(fourier_coeffsworld));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsworld(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fworldfilt4pt(k+1,1:flen) = ffiltwrap;
end;
%% recon filtered
fwfiltrecon4pt = SPHARMreconscalarcomplex(Fworldfilt4pt,L,sigma,theta,varphi,'normal');
if (norm(imag(fwfiltrecon4pt))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
Fsmoothr=real(fwfiltrecon4pt);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(3); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('filtered along 4pt L=%d',ll));
%% butterfly filter
% f(x,y) = x*exp(-(x^2+y^2)/2);
% stereographic: r = 2tan(theta), planar_theta = spherical_varphi;
hbflyx = zeros(size(theta));
hbflyy = zeros(size(theta));
butterfly_sigma = 0.01;
for r=1:size(theta,1)
    for c=1:size(theta,2)
        % stereogrphic projection
        rval = 2*tan(theta(r,c)/2);
        xval = rval*cos(varphi(r,c));
        yval = rval*sin(varphi(r,c));
        hbflyx(r,c) = xval*exp(-rval^2/butterfly_sigma);
    end;
end;
%% display butterfly
mF = min(hbflyx(:));
MF = max(hbflyx(:));
Fs = hbflyx - mF+1.0;  % now a function value of zero means it lies on the sphere
%Fs = 1 + (Fe-mF)/(MF-mF);
figure
trisurf(tri,x.*Fs,y.*Fs,z.*Fs,Fs);  % plot F(P)*P in Foley's notation
%title('Butterfly filter, looking down from north pole');
axis equal
axis off;
daspect([1 1 1]);
view(0,90); % down from north pole
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%% FIR approximation to butterfly
% create a uniform sampling grid 
% alpha is rotation around z axis, i.e., along lines of latitude
% beta is rotation around new x axis
% gamma is rotation around new z axis
nsample = 1;
betavals = []; alphvals = []; gammavals = [];
betavals(nsample) = 0; alphavals(nsample) = 0; gammavals(nsample) = 0;
hvalsx = []; hvalsy = [];
hvalsx(nsample) = 2*tan(0/2)*cos(0)*exp(-2*tan(0/2)/butterfly_sigma);
hvalsy(nsample) = 2*tan(0/2)*sin(0)*exp(-2*tan(0/2)/butterfly_sigma);
Ngrid = 11;  % 3 circles
theta_span = pi/8;
for g = 1:Ngrid; 
    thetaval = g/Ngrid * theta_span;
    Nsamples_per_circ = 2*g+1; % equiangular
    for v = 0:Nsamples_per_circ-1
       varphival = v/Nsamples_per_circ*2*pi;
       % convert from theta,varphi to alpha,beta,gamma
       nsample = nsample+1;
       alphavals(nsample)=varphival; betavals(nsample)=thetaval; gammavals(nsample)=-varphival;
       hvalsx(nsample) = 2*tan(thetaval/2)*cos(varphival)*exp(-2*tan(thetaval/2)/butterfly_sigma);
       hvalsy(nsample) = 2*tan(thetaval/2)*sin(varphival)*exp(-2*tan(thetaval/2)/butterfly_sigma);
    end;
end;
hvalsx=hvalsx/sum(abs(hvalsx));
hvalsy=hvalsy/sum(abs(hvalsy));
%% create sample points on plane from those on sphere
northpole = [0 0 1];
ptsonxy = zeros(nsample,2);
xvecloc = ptsonxy;
yvecloc = xvecloc;
for k=1:nsample
    R = euleranglestorotmatrix(alphavals(k),betavals(k),gammavals(k));
    ptvec = northpole*R;
    ptsonxy(k,:)=ptvec(1:2);
    %xvecproj = [1 0 0]*R;
    %xvecloc(k,:) = xvecproj(1:2);
    xvecloc(k,:) = -1*R(1,1:2);  % because we are looking down
    %yvecproj = [0 1 0]*R;
    %yvecloc(k,:) = yvecproj(1:2);
    yvecloc(k,:) = R(2,1:2);
end;
figure,
plot(ptsonxy(:,1),ptsonxy(:,2),'rx'),
xlabel('x'),ylabel('y');
axis equal
%axis off;
%hold on;
%quiver(ptsonxy(:,1),ptsonxy(:,2),0.05*xvecloc(:,1),0.05*xvecloc(:,2),0);
%quiver(ptsonxy(:,1),ptsonxy(:,2),0.05*yvecloc(:,1),0.05*yvecloc(:,2),0);
%hold off;
%title('X-Y axes projected on to the xy plane');
%% create fir filter
charok = 1;
L = BW-1;
for k=0:L-1
    firfiltx{k+1} = hvalsx(1)*eye(2*k+1);
    firfilty{k+1} = hvalsy(1)*eye(2*k+1);
end;
for n = 1:nsample
    DC = dmatrixbeta(L,betavals(n));
    if (n>1)
        for k = 0:L-1
            ind = -k:k;
            D = exp(j*ind'*alphavals(n))*exp(j*ind*gammavals(n)) .* DC{k+1};
            t = real(trace(D));
            
            % rotation angle around axis comes from the trace of
            % rep for k = 1
            if (k>1)
                u = sin(phi/2);
                if (u==0)
                    fprintf(1,'n=%d, k=%d phi=%f\n',n,k,phi);
                end;
                ch = sin((k+0.5)*phi)/sin(phi/2);
            else
                if k==1
                    phi = acos((t-1)/2);
                    if (phi == 0)
                        ch = 2*k+1;
                    else
                        ch = sin((k+0.5)*phi)/sin(phi/2);
                    end;
                else
                    ch = 1;
                end;
            end;
            if (abs (t-ch) > 1e-6)
                fprintf(1,'Error: rep for %d does not have right character\n',k);
                charok = 0;
            end;
            firfiltx{k+1} = firfiltx{k+1} + hvalsx(n)*D;
            firfilty{k+1} = firfilty{k+1} + hvalsy(n)*D;
        end;
    end;
end;
%% impulse response
Ffirfiltx = zeros(size(fourier_coeffsimpulse));
Ffirfilty = zeros(size(fourier_coeffsimpulse));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsimpulse(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*firfiltx{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Ffirfiltx(k+1,1:flen) = ffiltwrap;
    ffiltvec = funwrap*firfilty{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Ffirfilty(k+1,1:flen) = ffiltwrap;
end;
%% recon filtered x
ffiltreconx = SPHARMreconscalarcomplex(Ffirfiltx,L,sigma,theta,varphi,'normal');
if (norm(imag(ffiltreconx))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
%subplot(1,2,1);
Fsmoothr=real(ffiltreconx);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight headlight;
colormap gray
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('FIR approximation to butterfly in east-west direction L=%d',ll));
%% recon filtered y
ffiltrecony = SPHARMreconscalarcomplex(Ffirfilty,L,sigma,theta,varphi,'normal');
if (norm(imag(ffiltrecony))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
%subplot(1,2,2);
Fsmoothr=real(ffiltrecony);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
title(sprintf('FIR approximation to butterfly in north-south direction L=%d',ll));
%% filter the world with this approximation
Fwfirfiltx = zeros(size(fourier_coeffsworld));
Fwfirfilty = zeros(size(fourier_coeffsworld));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsworld(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*firfiltx{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fwfirfiltx(k+1,1:flen) = ffiltwrap;
    ffiltvec = funwrap*firfilty{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fwfirfilty(k+1,1:flen) = ffiltwrap;
end;
%% recon filtered world map in X direction
fwfiltreconx = SPHARMreconscalarcomplex(Fwfirfiltx,L,sigma,theta,varphi,'normal');
if (norm(imag(fwfiltreconx))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
%subplot(1,3,1);
Fsmoothr=abs(real(fwfiltreconx));
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(-165,30); % default 3D view;
%camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('world, filtered 0 deg longitude orientation'));
%% recon filtered world map in the y direction
fwfiltrecony = SPHARMreconscalarcomplex(Fwfirfilty,L,sigma,theta,varphi,'normal');
if (norm(imag(fwfiltrecony))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure,
%subplot(1,3,2);
Fsmoothr=abs(real(fwfiltrecony));
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(-165,30); % default 3D view;
%camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('world, filtered 90 degree longitude'));
figure; %subplot(1,3,3);
Fsmoothr=abs(real(fwfiltrecony)) + abs(real(fwfiltreconx));
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(-165,30); % default 3D view;
%camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('world, sum of both filters in |:|'));
%% dilate x filter in angle
dilation=2;  % 2x angles
charok = 1;
L = BW-1;
for k=0:L-1
    firfiltxd{k+1} = hvalsx(1)*eye(2*k+1);
    firfiltyd{k+1} = hvalsy(1)*eye(2*k+1);
end;
for n = 1:nsample
    DC = dmatrixbeta(L,dilation*betavals(n));
    if (n>1)
        for k = 0:L-1
            ind = -k:k;
            D = exp(j*ind'*dilation*alphavals(n))*exp(j*ind*dilation*gammavals(n)) .* DC{k+1};
            t = real(trace(D));
            
            % rotation angle around axis comes from the trace of
            % rep for k = 1
            if (k>1)
                u = sin(phi/2);
                if (u==0)
                    fprintf(1,'n=%d, k=%d phi=%f\n',n,k,phi);
                end;
                ch = sin((k+0.5)*phi)/sin(phi/2);
            else
                if k==1
                    phi = acos((t-1)/2);
                    if (phi == 0)
                        ch = 2*k+1;
                    else
                        ch = sin((k+0.5)*phi)/sin(phi/2);
                    end;
                else
                    ch = 1;
                end;
            end;
            if (abs (t-ch) > 1e-6)
                fprintf(1,'Error: rep for %d does not have right character\n',k);
                charok = 0;
            end;
            firfiltxd{k+1} = firfiltxd{k+1} + hvalsx(n)*D;
            firfiltyd{k+1} = firfiltyd{k+1} + hvalsy(n)*D;
        end;
    end;
end;
%% impulse response
Ffirfiltxd = zeros(size(fourier_coeffsimpulse));
Ffirfiltyd = zeros(size(fourier_coeffsimpulse));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsimpulse(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*firfiltxd{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Ffirfiltxd(k+1,1:flen) = ffiltwrap;
    ffiltvec = funwrap*firfiltyd{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Ffirfiltyd(k+1,1:flen) = ffiltwrap;
end;
%% recon filtered inpuluse response in x
ffiltreconxd = SPHARMreconscalarcomplex(Ffirfiltxd,L,sigma,theta,varphi,'normal');
if (norm(imag(ffiltreconxd))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
%subplot(1,2,1);
Fsmoothr=real(ffiltreconxd);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('FIR approximation to butterfly in east-west direction L=%d',ll));
%% recon filtered impulse response on y direction
ffiltreconyd = SPHARMreconscalarcomplex(Ffirfiltyd,L,sigma,theta,varphi,'normal');
if (norm(imag(ffiltreconyd))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
%subplot(1,2,2);
Fsmoothr=real(ffiltreconyd);
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('FIR approximation to butterfly in north-south direction L=%d',ll));
%% filter the world with this approximation
Fwfirfiltxd = zeros(size(fourier_coeffsworld));
Fwfirfiltyd = zeros(size(fourier_coeffsworld));
for k=0:L-1
    flen = 2*k+1;
    fvec = fourier_coeffsworld(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*firfiltxd{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fwfirfiltxd(k+1,1:flen) = ffiltwrap;
    ffiltvec = funwrap*firfiltyd{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fwfirfiltyd(k+1,1:flen) = ffiltwrap;
end;
%% recon filtered world map after dilation
fwfiltreconxd = SPHARMreconscalarcomplex(Fwfirfiltxd,L,sigma,theta,varphi,'normal');
if (norm(imag(fwfiltreconxd))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
%subplot(1,3,1);
Fsmoothr=abs(real(fwfiltreconxd));
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(-165,30); % default 3D view;
%camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('world, filtered 0 deg longitude orientation dilated 2x'));
%% recon filtered world map in y after dilation
fwfiltreconyd = SPHARMreconscalarcomplex(Fwfirfiltyd,L,sigma,theta,varphi,'normal');
if (norm(imag(fwfiltreconyd))>1e-6)
    disp('Error! Imaginary part is too large');
    return;
end;
figure;
%subplot(1,3,2);
Fsmoothr=abs(real(fwfiltreconyd));
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(-165,30); % default 3D view;
%camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('world, filtered 90 degree longitude dilated 2x'));
%subplot(1,3,3);
figure;
Fsmoothr=abs(real(fwfiltreconyd)) + abs(real(fwfiltreconxd));
mF = min(Fsmoothr(:));
MF = max(Fsmoothr(:));
Fsf = 1 + (Fsmoothr-mF);
trisurf(tri,x.*Fsf,y.*Fsf,z.*Fsf,Fsf);
axis equal;
axis off;
daspect([1 1 1]);
view(-165,30); % default 3D view;
%camlight headlight;
%axis vis3d off;
lighting phong;
material dull; %shiny;
shading interp;
%title(sprintf('world, sum of both filters in |:| dilated 2x'));
%% Download  a sample outer cortical surface outersurface.obj and read it with mni_getmesh.m.
[tri,coord,nbr,normal]=mni_getmesh('outersurface.obj');
figure,trisurf(tri,coord(1,:),coord(2,:),coord(3,:)); %,title('data');
axis tight; 
shading interp; 
colormap bone;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight; % light placed at camera position
lighting gouraud; % match what Chung used
%alpha(0.9); % what Chung has
%%
%Run SPHARMsmoothcomplex.m to get the complex weighted-SPHARM representation. 
%this uses 80th degree weighted-SPHARM with bandwidth 0.0001 

L=63;
sigma=0.0001;
[b_coord_smooth_c,b_fourier_coeff_c]= SPHARMsmoothcomplex(coord,L,sigma);
figure,trisurf(tri,b_coord_smooth_c(1,:),b_coord_smooth_c(2,:),b_coord_smooth_c(3,:));
%title('complex smooth reconstruction');
axis tight; 
shading interp; 
colormap bone;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight; % light placed at camera position
lighting gouraud; % match what Chung used
%alpha(0.9); % what Chung has
%% now show reconstruction for complex-valued coeffs
load unitsphere.mat;
sphere.vertices=coord';
sphere.faces=tri; %+1;
L=63;
sigma=0.0001;
resurf= SPHARMrepresentcomplex(sphere,b_fourier_coeff_c,L);
cr = resurf.vertices'; 
%tr = resurf.faces;
figure,trisurf(tri,cr(:,1),cr(:,2),cr(:,3));
title('complex reconstruction');
axis tight; 
shading interp; 
colormap bone;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight; % light placed at camera position
lighting gouraud; % match what Chung used
%alpha(0.9); % what Chung has
%% 4pt filtering of filter
Fbfilt4pt.x = zeros(size(b_fourier_coeff_c.x));
Fbfilt4pt.y = zeros(size(b_fourier_coeff_c.y));
Fbfilt4pt.z = zeros(size(b_fourier_coeff_c.z));
for k=0:L-1
    flen = 2*k+1;
    fvec = b_fourier_coeff_c.x(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfilt4pt.x(k+1,1:flen) = ffiltwrap;
    fvec = b_fourier_coeff_c.y(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfilt4pt.y(k+1,1:flen) = ffiltwrap;
    fvec = b_fourier_coeff_c.z(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfilt4pt.z(k+1,1:flen) = ffiltwrap;
end;
%% now show reconstruction for complex-valued coeffs
load unitsphere.mat;
sphere.vertices=coord';
sphere.faces=tri; %+1;
L=63;
sigma=0.0001;
resurf= SPHARMrepresentcomplex(sphere,Fbfilt4pt,L);
cr = resurf.vertices'; 
%tr = resurf.faces;
figure,trisurf(tri,cr(:,1),cr(:,2),cr(:,3));
%title('complex filtered reconstruction');
axis tight; 
shading interp; 
colormap bone;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight; % light placed at camera position
lighting gouraud; % match what Chung used
%alpha(0.9); % what Chung has
%% repeat 4pt filtering of filter
Fbfilt4pt2.x = zeros(size(b_fourier_coeff_c.x));
Fbfilt4pt2.y = zeros(size(b_fourier_coeff_c.y));
Fbfilt4pt2.z = zeros(size(b_fourier_coeff_c.z));
for k=0:L-1
    flen = 2*k+1;
    fvec = b_fourier_coeff_c.x(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1}*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfilt4pt2.x(k+1,1:flen) = ffiltwrap;
    fvec = b_fourier_coeff_c.y(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1}*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfilt4pt2.y(k+1,1:flen) = ffiltwrap;
    fvec = b_fourier_coeff_c.z(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmat4pt{k+1}*Dfiltmat4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfilt4pt2.z(k+1,1:flen) = ffiltwrap;
end;
%% now show reconstruction for complex-valued coeffs
load unitsphere.mat;
sphere.vertices=coord';
sphere.faces=tri; %+1;
L=63;
sigma=0.0001;
resurf= SPHARMrepresentcomplex(sphere,Fbfilt4pt2,L);
cr = resurf.vertices'; 
%tr = resurf.faces;
figure,trisurf(tri,cr(:,1),cr(:,2),cr(:,3));
%title('complex filtered reconstruction');
axis tight; 
shading interp; 
colormap bone;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight; % light placed at camera position
lighting gouraud; % match what Chung used
%alpha(0.9); % what Chung has
%% create a simple 5 pt unsharp masking filter both lat and long lines filter
% alpha is rotation around z axis, i.e., along lines of latitude
% beta is rotation around new x axis
% gamma is rotation around new z axis
beta1 = pi/32; alpha1=pi/2; gamma1=-pi/2;  % latitude only
beta2 = pi/32; alpha2 = 0; gamma2=0; % longitude only
%% check if the reps have the right character
DC1 = dmatrixbeta(L,beta1);  % identify
DC2 = dmatrixbeta(L,beta2);
charok = 1;
for k=0:L-1
    Dfiltmatusm4pt{k+1} = eye(2*k+1);
end;
for k=0:L-1
    ind = -k:k;
    D1 = exp(j*ind'*alpha1)*exp(j*ind*gamma1) .* DC1{k+1};
    D2 = exp(j*ind'*alpha2)*exp(j*ind*gamma2) .* DC2{k+1};
    
    t1 = real(trace(D1));
    
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)
        ch = sin((k+0.5)*phi)/sin(phi/2);
    else
        if k==1
            phi = acos((t1-1)/2);
            if (phi == 0)
                ch = 2*k+1;
            else
                ch = sin((k+0.5)*phi)/sin(phi/2);
            end;
        else
            ch = 1;
        end;
    end;
    if (abs (t1-ch) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
    t2 = real(trace(D2));
    
    % rotation angle around axis comes from the trace of
    % rep for k = 1
    if (k>1)
        ch2 = sin((k+0.5)*phi2)/sin(phi2/2);
    else
        if k==1
            phi2 = acos((t2-1)/2);
            if (phi2 == 0)
                ch2 = 2*k+1;
            else
                ch2 = sin((k+0.5)*phi2)/sin(phi2/2);
            end;
        else
            ch2 = 1;
        end;
    end;
    if (abs (t2-ch2) > 1e-6)
        fprintf(1,'Error: rep for %d does not have right character\n',k);
        charok = 0;
    end;
    Dfiltmatusm4pt{k+1}=Dfiltmatusm4pt{k+1}-0.65*(4*eye(size(D1))+D1+D1'+D2+D2')/8;
end;
if (charok)
    fprintf(1,'Character test passed\n');
end;
%% 4pt filtering of filter
Fbfiltusm4pt.x = zeros(size(b_fourier_coeff_c.x));
Fbfiltusm4pt.y = zeros(size(b_fourier_coeff_c.y));
Fbfiltusm4pt.z = zeros(size(b_fourier_coeff_c.z));
for k=0:L-1
    flen = 2*k+1;
    fvec = b_fourier_coeff_c.x(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmatusm4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfiltusm4pt.x(k+1,1:flen) = ffiltwrap;
    fvec = b_fourier_coeff_c.y(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmatusm4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfiltusm4pt.y(k+1,1:flen) = ffiltwrap;
    fvec = b_fourier_coeff_c.z(k+1,1:flen);
    funwrap = [fliplr(fvec(k+2:end)),fvec(1),fvec(2:k+1)];
    ffiltvec = funwrap*Dfiltmatusm4pt{k+1};
    ffiltwrap = [ffiltvec(k+1),ffiltvec(k+2:end),fliplr(ffiltvec(1:k))]; %imyfftshift(frotvec);
    Fbfiltusm4pt.z(k+1,1:flen) = ffiltwrap;
end;
%% now show reconstruction for complex-valued coeffs
load unitsphere.mat;
sphere.vertices=coord';
sphere.faces=tri; %+1;
L=40;
sigma=0.0001;
resurf= SPHARMrepresentcomplex(sphere,Fbfiltusm4pt,L);
cr = resurf.vertices'; 
mc = min(cr(:));
MC = max(cr(:));
cr2 = 1+(cr-mc)/(MC-mc);
%tr = resurf.faces;
figure,trisurf(tri,cr(:,1),cr(:,2),cr(:,3));
title('complex unshap filtered reconstruction');
axis tight; 
shading interp; 
colormap bone;
axis off;
daspect([1 1 1]);
view(0,90); % default 3D view;
camlight; % light placed at camera position
lighting gouraud; % match what Chung used
%alpha(0.9); % what Chung has