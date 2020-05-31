%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate what the field looks like at the detection plane 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code generates an input field |U> and computes the optical overlap
% between it and the basis mode \Phi_n. 
% Can be used for diagnostic purposes to check if the simulated intensity
% matches what is observed with a camera.
% Adjust the parameters below as required

% Make square coordinate system (in units of mm)
H = 1000; % number of points
dx = 8e-3; % pixel size
x = dx.*(-H/2:(H/2-1)); % x coordinates (same as y)
[X,Y] = meshgrid(x,-x); % matrix Cartesian coordinates
[Phi,R] = cart2pol(X,Y); % matrix Polar coordinates

% Parameters of optical system
f = 1000; % focal length of Fourier lens
lambda = 633e-6; % wavelength of light

% Make input field
w0 = 0.5; % beam size
% LG beam
L = [5,-5]; P = [3,3]; % LG mode indices
c = [1,1]; % expansion coefficients (see LG.m)
Mode_in = LG(R,Phi,P,L,c,w0);
% HG beam
% N = [4,2]; M = [2,4]; % HG mode indices
% c = [1,1]; % expansion coefficients (see HG.m)
% Mode_in = HG(X,Y,N,M,c,w0);

% Make basis element (\Phi_n)
w00 = 0.5;
% LG basis
LL = [5]; PP = [3]; % LG mode indices
BasisEl = LG(R,Phi,PP,LL,1,w00);
% HG basis
% NN = 4; MM = 2; % HG mode indices
% BasisEl = HG(X,Y,NN,MM,1,w00);

% Make discrete optical Fourier transform matrix
L = H*dx; % physical side length of image
Lk = lambda*f/dx; %physical side length at Fourier plane 
dk = lambda*f/L; %sample size at Fourier plane
k = -Lk/2:dk:Lk/2-dk; %spatial frequency coordinate system
k = k./dk^2; % rescale for optical spatial frequencies
DFT = exp(-1i*2*pi/H).^(k'*x); % Fourier transform matrix

% Field at detection plane is the Fourier transform of the input field
% modulated by the basis mode
FT_field = DFT'*(Mode_in.*conj(BasisEl)).'*DFT;

% Plots: an on-axis bright spot indicates the basis element forms part of
% the input field. An on-axis null indicates it does not form part of the
% input field
Q = 300; % image size
figure('color','w','units','points','position',[50 100 3*Q Q]);

subplot(1,3,1);
imagesc(abs(Mode_in).^2);
text(H/40,H/20,'Input field |U\rangle','FontSize',20);
set(gca,'units','points','position',[0 0 Q Q],'visible','off');

subplot(1,3,2);
imagesc(abs(BasisEl).^2);
text(H/40,H/15,'Basis element \langle\Phi_n|','FontSize',20);
set(gca,'units','points','position',[Q 0 Q Q],'visible','off');

subplot(1,3,3);
imagesc(abs(FT_field).^2);
text(H/40,H/15,'At detector \langle\Phi_n|U\rangle','FontSize',20);
set(gca,'units','points','position',[2*Q 0 Q Q],'visible','off');

% Phase plots as inserts
axes('pos',[.23 .7 .1 .3])
imagesc(angle(Mode_in)); xticks([]); yticks([]);
axes('pos',[.56 .7 .1 .3])
imagesc(angle(BasisEl)); xticks([]); yticks([]);
axes('pos',[.9 .7 .1 .3])
imagesc(angle(FT_field)); xticks([]); yticks([]);
