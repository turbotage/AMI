% H-scan std map and contrast computations
% ------------------------------------------------------------------------
% Pre Filters (choose dataset)
%-----------------------------
% Dataset 2
RFMATlo=squeeze(Bmodesrgblo(:,:,3,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,1,:));

RF_MATlo_f = mybandpass(RFMATlo,4,[3,30],500);
RF_MATlo_f = mymedfilt(RF_MATlo_f, [10,3]);

RFMAThi_f = mybandpass(RF_MAThi,4,[3,30],500);
RFMAThi_f = mymedfilt(RFMAThi_f, [10,3]);

%% ----------------------------
%  Calculate contrast in images
% -----------------------------
stdratiolo = roi_metric(RF_MATlo_f, roimat);
stdratiohi = roi_metric(RFMAThi_f, roimat);

fprintf('\nContrast for low pass: %.3f \n',stdratiolo)
fprintf('Contrast for high pass: %.3f \n',stdratiohi)

%% Dataset 3
RF_MATlo=squeeze(Bmodesrgblo(:,:,3,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,1,:));

RF_MATlo_f = mybandpass(RF_MATlo,4,[5,25],1000);
RF_MATlo_f = mymedfilt(RF_MATlo_f, [10,3]);

RF_MAThi_f = mybandpass(RF_MAThi,4,[5,25],1000);
RF_MAThi_f = mymedfilt(RF_MAThi_f, [10,3]);

disp('Pre std mask filters of H-scan images')

%% ----------------------------
% Create and plot std masks
%------------------------------
perc=[1,99];
draw_std2(RF_MATlo_f,RF_MAThi_f,perc);
