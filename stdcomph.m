% H-scan std map and contrast computations
% ------------------------------------------------------------------------
% Pre Filters (choose dataset)
%-----------------------------
% Dataset 2
RFMATlo=squeeze(Bmodesrgblo(:,:,3,:));
RFMAThi=squeeze(Bmodesrgbhi(:,:,1,:));

RFMATlo_f = mybandpass(RFMATlo,4,[3,30],500);
RFMATlo_f = mymedfilt(RFMATlo_f, [10,3]);

RFMAThi_f = mybandpass(RFMAThi,4,[3,30],500);
RFMAThi_f = mymedfilt(RFMAThi_f, [10,3]);

%% ----------------------------
%  Calculate contrast in images
% -----------------------------
stdratiolo = roi_metric(RFMATlo_f, roimat);
stdratiohi = roi_metric(RFMAThi_f, roimat);

fprintf('\nContrast for low pass: %.3d \n',stdratiolo)
fprintf('Contrast for high pass: %.3d \n',stdratiohi)

%% Dataset 3
RFMATlo=squeeze(Bmodesrgblo(:,:,3,:));
RFMAThi=squeeze(Bmodesrgbhi(:,:,1,:));

RFMATlo_f = mybandpass(RFMATlo,4,[5,25],1000);
RFMATlo_f = mymedfilt(RFMATlo_f, [10,3]);

RFMAThi_f = mybandpass(RFMAThi,4,[5,25],1000);
RFMAThi_f = mymedfilt(RFMAThi_f, [10,3]);

disp('Pre std mask filters of H-scan images')

%% ----------------------------
% Create and plot std masks
%------------------------------
perc=[1,99];
draw_std2(RFMATlo_f,RFMAThi_f,perc);


