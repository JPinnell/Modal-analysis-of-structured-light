%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate a modal decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code generates a structured light field |U> and performs a modal
% decomposition in the basis Phi_n (taken to be the LG basis but you can 
% change to another basis if desired).
% The example below is that given in our paper "Modal analysis of
% structured light: a practical tutorial"
% Uses the in-built fast fourier transform (fft2) to compute the overlap

% Make square coordinate system (in units of mm)
H = 1000; % number of points
dx = 8e-3; % pixel size
x = dx.*(-H/2:(H/2-1)); % x coordinates (same as y)
[X,Y] = meshgrid(x,-x); % matrix Cartesian coordinates
[Phi,R] = cart2pol(X,Y); % matrix Polar coordinates

% Make the initial field (superposition of LG modes in Eq. 29)
w0 = 1; % beam size of the field
P = [0,1,2,0,1,2]; % radial indices
L = [1,-1,1,0,1,-1]; % azimuthal indices
c = [0.5 0.25 1 0.5 0.25 1].*exp(1i.*[pi pi/4 -pi/2 0 pi/4 -pi/2]); % expansion coefficients
InputField = LG(R,Phi,P,L,c,w0); % see LG.m
ReconField = zeros(H,H); % initialise reconstructed field

% Do the decomposition for a predetermined set of modes
PMax = 2; % maximum decomp for radial index
LMax = 1; % maximum decomp for azimuthal index
c_meas = zeros(PMax+1,2*LMax+1); % initialise
w00 = w0; % beam size of basis modes
for l = -LMax:LMax
    for p = 0:PMax
        BasisEl = LG(R,Phi,p,l,1,w00);
        overlap = fftshift(fft2(InputField.*conj(BasisEl))); % Fourier transform modulated field
        c_meas(l+LMax+1,p+1) = overlap(H/2+1,H/2+1); % extract on-axis portion
        ReconField = ReconField + c_meas(l+LMax+1,p+1).*BasisEl; % recondstructed field
    end
end
c_meas = c_meas./norm(c_meas); % manually normalise measured expansion coefficients

% Plots
Q = 200;
figure('color','w','units','points','position',[50 50 2*Q 2*Q]);

subplot(2,2,1);
imagesc(abs(InputField).^2);
text(H/30,H/20,'Input','FontSize',14);
set(gca,'units','points','position',[0 Q Q Q],'visible','off');

subplot(2,2,2);
imagesc(abs(ReconField).^2);
text(H/30,H/20,'Reconstructed','FontSize',14);
set(gca,'units','points','position',[Q Q Q Q],'visible','off');

subplot(2,2,3);
bar(abs(c_meas(:))); title('\rho_n');
set(gca,'units','points','position',[Q/10 Q/10 0.8*Q 0.8*Q]);

subplot(2,2,4);
bar(angle(c_meas(:))); title('\phi_n');
set(gca,'units','points','position',[Q+Q/10 Q/10 0.8*Q 0.8*Q]);

% Phase plots as inserts
axes('pos',[.35 .85 .15 .15])
imagesc(angle(InputField)); xticks([]); yticks([]);
axes('pos',[.85 .85 .15 .15])
imagesc(angle(ReconField)); xticks([]); yticks([]);