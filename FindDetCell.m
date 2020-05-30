function DetCellMask = FindDetCell(Img,PixelSize,dmin,threshold)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function help
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DetCellMask] = FindDetCell(Img,PixelSize,dmin)
% v1 J.Pinnell 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Descrition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the detection cell mask corresponding to the
% on-axis region when using a camera for optical modal decomposition. It
% does this by fitting a Gaussian to the beam and finding the coordinates
% of the maximum.
% NB: The on-axis intensity is then the average of the intensity in this 
% region, use -> mean(Img.*DetCellMask)
% Details found in Pinnell et. al., "Modal analysis of structured light: a practical tutorial"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Img - camera image of the focused beam at the Fourier/detection plane
% PixelSize - pixel size of the camera (assumes square pixels)
% dmin - minimum resolvable distance of the optics (angular resolution)
% threshold - thresholding the Gaussian fit (~0.3 often works well)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DetCellMask - the binary mask corresponding to the detection cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H = 1000;
% PixelSize = 8e-3; x = PixelSize.*(-H/2:H/2-1); [X,Y] = meshgrid(x,-x); 
% w0 = 0.5; dmin = w0; threshold = 0.3; shift = randn(2,1);
% Img = exp(-((X-shift(1)).^2 + (Y-shift(2)).^2)./w0.^2);
% DetCellMask = FindDetCell(Img,PixelSize,dmin,threshold);

[V,H] = size(Img);

avgX = mean(Img,1); avgX = avgX./max(avgX); % average image in x direction
avgY = mean(Img,2)'; avgY = avgY./max(avgY); % average image in y direction

x = 1:H; % Camera column (x) dimensions
y = 1:V; % Camera row (y) dimensions

polyX = polyfit(x(avgX > threshold),log(avgX(avgX > threshold)),2); % x average fit
polyY = polyfit(y(avgY > threshold),log(avgY(avgY > threshold)),2); % y average fit

fitX = exp(polyval(polyX,x)); % fitted x curve
fitY = exp(polyval(polyY,y)); % fitted y curve

[~,Cx] = max(fitX); % maximum col (x) of fitted curve
[~,Cy] = max(fitY); % maximum row (y) of fitted curve

[XX,YY] = meshgrid(x,y); 
DetCellMask = (((XX-Cx).^2+(YY-Cy).^2)<(dmin/2/PixelSize)^2);

% Plots (comment out if you don't want to see it everytime)
Q = 200;
figure('color','w','units','points','position',[50 50 400 400]);

subplot(2,2,1);
imagesc(DetCellMask);
text(H/30,V/10,'Mask','FontSize',20);
set(gca,'units','points','position',[0 Q Q Q],'visible','off');

subplot(2,2,2);
hold on;
plot(x,avgX,'-r','LineWidth',4);
plot(x,fitX,'-k','LineWidth',2);
line([Cx,Cx],[0,max(fitX)],'Color',[1,0,0],'LineWidth',2,'LineStyle','--');
xlim([1,H]); ylim([0,1]);
hold off;
text(H/30,0.9,'X fit','FontSize',20);
set(gca,'units','points','position',[Q Q Q Q],'visible','off');

subplot(2,2,3);
hold on;
plot(-avgY,-y,'-b','LineWidth',4);
plot(-fitY,-y,'-k','LineWidth',2);
line([0,-max(fitY)],[-Cy,-Cy],'Color',[1,0,0],'LineWidth',2,'LineStyle','--');
xlim([-1,0]); ylim([-V,-1]);
hold off;
text(-0.96,-V/10,'Y fit','FontSize',20);
set(gca,'units','points','position',[0 0 Q Q],'visible','off');

subplot(2,2,4);
imagesc(Img);
viscircles([Cx Cy], dmin/2/PixelSize, 'Color',[1,0,0],'LineWidth',2); % might not work if you dont have image processing toolbox
text(H/30,V/10,'Image','FontSize',20);
set(gca,'units','points','position',[Q 0 Q Q],'visible','off');

end