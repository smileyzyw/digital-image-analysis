clear;
close all;
clc;

% Img=imread('three.bmp');      % example works well
% Img=imread('twoCells.bmp');     % example works well
Img=imread('vessel.bmp');     % example does NOT work well
U=Img(:,:,1);
I=double(U);

% get the size
[nrow,ncol] =size(U);
ic=nrow/2;
jc=ncol/2;
r=30;
phi_0 = sdf2circle(nrow,ncol,ic,jc,r);

delta_t = 10;   %   time step
mu = 0.001;      %   distRictTerm coefficient
nu = 1;         %   areaTerm coefficient
lambda=3;       %   legthTerm coefficient

epsilon=1.5;    %   used in Dirac function 

% iteration should begin from here
phi=phi_0;
numIter = 1;
for k=1:200
    phi = evolution(I, phi, mu, nu, lambda, delta_t, epsilon, numIter);  
    if mod(k,10) == 0
        pause(0.05);
        figure(2);
        imagesc(uint8(I));colormap(gray)
        hold on;
        plotLevelSet(phi,0,'r');
    end    
end

