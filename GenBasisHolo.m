function [Holo] = GenBasisHolo(Mode,X,Y,Nx,Ny)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Holo] = GenBasisHolo(Mode,X,Y,Nx,Ny)
% v1 J.Pinnell 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Descrition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function makes the phase-only hologram corresponding to a given
% basis mode. It uses exact complex amplitude modulation as outlined in
% Bolduc, Optics Letters 38(18) (2003).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mode - transverse (complex) electric field profile of the basis function
% X,Y - 2D coordinate system of the mode
% Nx,Ny - number of diffraction grating lines in x,y directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Holo - hologram to be displayed on the phase-only SLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = 1000; PixelSize = 8e-3; x = PixelSize.*(-H/2:H/2-1); 
% [X,Y] = meshgrid(x,-x); [Phi,R] = cart2pol(X,Y); Nx = 50; Ny = 50;
% Mode = LG(R,Phi,[5,5],[5,-5],[1,1],0.9); % see LG.m
% imagesc(GenBasisHolo(Mode,X,Y,Nx,Ny)); colormap gray;

[V,H] = size(Mode);
dx = abs(X(1,2)-X(1,1)); % pixel size in x direction
dy = abs(Y(2,1)-Y(1,1)); % pixel size in x direction
Gx = Nx/(H*dx); Gy = Ny/(V*dy); % spatial frequencies in units of lines per screen

A = abs(Mode); A = A./max(max(A)); % mode amplitude
Phase = angle(Mode); % phase of mode

load('SincInv.mat'); % load inverted sinc
aux = round(A.*(length(SincInv)-1)+1); % scale so that amplitudes map to amplitudes
F = 1 + SincInv(aux); % modified amplitude
Phi = mod(Phase-pi.*F+2*pi.*(Gx.*X+Gy.*Y),2*pi); % modified phase

Holo = F.*Phi; % hologram phase map 

Holo = Holo./max(Holo(:)); % Normalise
%Holo = uint8(Holo*255); % Scale pixels to [0,255] for 8 bit SLM screen
end