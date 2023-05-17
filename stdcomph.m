% AMI Project 2022 - Rebecca
% Script generationg standard deviation mask and contrast values
% - H-scan data from dataset 2 and 3 respectively
%----------------------------
% H-scan std map and contrast computations for both dataset 2 and 3
% respectively. Chose dataset to visualized by first running its
% "Pre-filter" section followed by the "Create plot and std masks" at the
% bottom of this file.

%% ---------------------------
% Pre-filters Dataset 2
%-----------------------------

RF_MATlo=squeeze(Bmodesrgblo(:,:,1,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,3,:));

% Comparative statement for selecting the proper filtering based on
% stimulation frequency of the previously chosen dataset.

if strcmp(chosen_dataset,'1555') || strcmp(chosen_dataset,'1311')
    RF_MATlo_f = mybandpass(RF_MATlo,4,[3,20],500);
	RF_MAThi_f = mybandpass(RF_MAThi,4,[3,20],500);

else
	RF_MATlo_f = mybandpass(RF_MATlo,4,[.5,10],500);
	RF_MAThi_f = mybandpass(RF_MAThi,4,[.5,10],500);
end

chosen_dataset=strcat('Dataset 2\_',chosen_dataset);
disp(chosen_dataset)

%% -------------------------------------------------
%  Calculate contrast in images (only for Dataset 2)
% --------------------------------------------------
stdratiolo = roi_metric(RF_MATlo, roimat);
stdratiohi = roi_metric(RF_MAThi, roimat);

stdratiolo_f = roi_metric(RF_MATlo_f, roimat);
stdratiohi_f = roi_metric(RF_MAThi_f, roimat);

fprintf('\nContrast for low pass wo. bp: %.3f \n',stdratiolo)
fprintf('Contrast for high pass wo. bp: %.3f \n',stdratiohi)

fprintf('\nContrast for low pass with bp: %.3f \n',stdratiolo_f)
fprintf('Contrast for high pass with bp: %.3f \n',stdratiohi_f)

%% ---------------------------
% Pre-filters Dataset 3
%-----------------------------
chosen_dataset='Dataset 3 RF\_MAT';
disp(chosen_dataset)

RF_MATlo=squeeze(Bmodesrgblo(:,:,1,:));
RF_MAThi=squeeze(Bmodesrgbhi(:,:,3,:));

RF_MATlo_f = mybandpass(RF_MATlo,4,[4,30],1000);
RF_MAThi_f = mybandpass(RF_MAThi,4,[4,30],1000);

%% ----------------------------
% Create and plot std masks
%------------------------------
% wo. post band pass filtering (BP)

perc=[1,99];
draw_std2(chosen_dataset,RF_MATlo,RF_MAThi,perc,roimat);

%% with post band pass filtering (BP)

chosen_dataset1=strcat(chosen_dataset,' with bp filtering');
perc=[1,99];
draw_std2(chosen_dataset1,RF_MATlo_f,RF_MAThi_f,perc,roimat);
