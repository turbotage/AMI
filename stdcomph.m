% AMI Project 2022 - Rebecca
% H-scan std map and contrast computations

% ----------------------------
% Pre-filters (choose dataset)
%-----------------------------

% Dataset 2
RF_MATlo=squeeze(Bmodesrgblo(:,:,1,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,3,:));

if strcmp(chosen_dataset,'1555') || strcmp(chosen_dataset,'1311')
    RF_MATlo_f = mybandpass(RF_MATlo,4,[3,20],500);
	RF_MAThi_f = mybandpass(RF_MAThi,4,[3,20],500);

else
	RF_MATlo_f = mybandpass(RF_MATlo,4,[.5,10],500);
	RF_MAThi_f = mybandpass(RF_MAThi,4,[.5,10],500);
end


disp('Pre std mask filters of H-scan images')

% ----------------------------
%  Calculate contrast in images
% -----------------------------
stdratiolo = roi_metric(RF_MATlo, roimat);
stdratiohi = roi_metric(RF_MAThi, roimat);

stdratiolo_f = roi_metric(RF_MATlo_f, roimat);
stdratiohi_f = roi_metric(RF_MAThi_f, roimat);

fprintf('\nContrast for low pass wo. bp: %.3f \n',stdratiolo)
fprintf('Contrast for high pass wo. bp: %.3f \n',stdratiohi)

fprintf('\nContrast for low pass with bp: %.3f \n',stdratiolo_f)
fprintf('Contrast for high pass with bp: %.3f \n',stdratiohi_f)

%% Dataset 3
RF_MATlo=squeeze(Bmodesrgblo(:,:,1,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,3,:));

RF_MATlo_f = mybandpass(RF_MATlo,4,[4,30],1000);
RF_MAThi_f = mybandpass(RF_MAThi,4,[4,30],1000);

disp('Pre std mask filters of H-scan images')

%% ----------------------------
% Create and plot std masks
%------------------------------
% wo bp
perc=[1,99];
draw_std2(RF_MATlo,RF_MAThi,perc);

%% with bp
perc=[1,99];
draw_std2(RF_MATlo_f,RF_MAThi_f,perc);
