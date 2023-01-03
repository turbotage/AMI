% AMI Project 2022 - Rebecca
% H-scan std map and contrast computations
% ------------------------------------------------------------------------
% Pre Filters (choose dataset)
%-----------------------------
% Dataset 2
RF_MATlo=squeeze(Bmodesrgblo(:,:,1,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,3,:));

RF_MATlo_f = mybandpass(RF_MATlo,4,[3,30],500);
RF_MAThi_f = mybandpass(RF_MAThi,4,[3,30],500);
disp('Pre std mask filters of H-scan images')

% ----------------------------
%  Calculate contrast in images
% -----------------------------
stdratiolo = roi_metric(RF_MATlo_f, roimat);
stdratiohi = roi_metric(RF_MAThi_f, roimat);

fprintf('\nContrast for low pass: %.3f \n',stdratiolo)
fprintf('Contrast for high pass: %.3f \n',stdratiohi)

%% Dataset 3
RF_MATlo=squeeze(Bmodesrgblo(:,:,1,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,3,:));

RF_MATlo_f = mybandpass(RF_MATlo,4,[4,30],1000);

RF_MAThi_f = mybandpass(RF_MAThi,4,[5,30],1000);
disp('Pre std mask filters of H-scan images')

%% ----------------------------
% Create and plot std masks
%------------------------------
perc=[1,99];
draw_std2(RF_MATlo_f,RF_MAThi_f,perc);
