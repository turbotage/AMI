%  AMI Project 2022

%% load dataset 1 only image contraction at full activation vs full rest
% should be conducted on one frame with only the H-scan for both full
% activation and at rest.

% --------
% AT REST
% --------
frame = 10; % Chose a frame, for dataset 1 one is sufficient for the report.
line = 64; % Image line = 64 out or 128
load('C:\Users\Rebecca Viklund\Desktop\AMI project\AMI\Dataset_1\CG_rest_2_omg2.mat');

% ----------------------
% Pre-processing of data
% ----------------------
RF = double(RF(1:2177,:,frame)); %1:2177 corresp to 4cm depth in this dataset..

% TGC - Only needed for dataset 1. TGC = "time gain compensation"
% It amplifies the signal proportional to the US depth since it weakens as
% it interacts with tissue.
TGC_VECT = linspace(1,10,size(RF,1))';
TGC_MAT = repmat(TGC_VECT,[1 size(RF,2)]);

RF = RF.*TGC_MAT;

% Here a "beamforming" is conducted. This is only needed for dataset 1.
% This is due to the collection method "plane-wave imaging" used which
% requires synthetic focusing which is provided by beamforming techniques.

BF = beamformDAQ(RF,128);
BF = BF - repmat(mean(BF,1),[size(BF,1) 1]);

% Use the beamformed version of RF for all computations below
RF = BF;

RF = single(RF(:,:));
Bmodes = single(sqrt(abs(hilbert(RF(:,:)))));
shape = size(Bmodes);
% ----------------------
% H-scan 1a
%-----------------------
Fs = 35; % MHz Samplingfreq of the US not the signal
T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);

% GH low pass
b1 = 0.13;
ordlo = 8;
Hlo = hermiteH(ordlo, t./b1);
GHlo = exp(-(t./(b1)).^2).*Hlo;
GHlo = GHlo./sum(GHlo(:));
[pxxlo,flo] = pwelch(GHlo, hamming(512));
f_VECTlo = linspace(0,Fs/2,length(pxxlo));
p_NORMlo = 0.5*pxxlo./max(pxxlo);

% GH high pass
b2 = 0.135;
ordhi = 32;
Hhi = hermiteH(32, t./b2); % order 32
GHhi = exp(-(t./(b2)).^2).*Hhi;
GHhi = GHhi./sum(GHhi(:));
[pxxhi,fhi] = pwelch(GHhi,hamming(512));
f_VECThi = linspace(0,Fs/2,length(pxxhi));
p_NORMhi = 0.5*pxxhi./max(pxxhi);

% Plot one spectrum of one imageline for comparison with the GH spectra
imLineRF = double(squeeze(RF(:,line))); 
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));
p_NORM = sqrt(pxx)./max(sqrt(pxx));

figure(1);
subplot(2,1,1)
plot(t,GHlo,'b', 'LineWidth', 1.5);
title(['GH_{',num2str(ordlo),'}'])
ylim([min(GHlo,[],'all')*1.25 max(GHlo,[],'all')*1.25]);
subplot(2,1,2)
plot(t,GHhi,'r', 'LineWidth', 1.5);
title(['GH_{',num2str(ordhi),'}'])
ylim([min(GHhi,[],'all')*1.25 max(GHhi,[],'all')*1.25]);
sgtitle('Gaussian weighted Hermite polynomials')

figure(2); clf; hold on;
plot(f_VECT, p_NORM,'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECTlo, p_NORMlo, 'b', 'LineWidth', 1.5);
plot(f_VECThi, p_NORMhi, 'r', 'LineWidth', 1.5);
sgtitle(['PSD for frame ',num2str(frame),' and line ',num2str(line)]);
xlabel('Frequency f [MHz]'); ylabel('Amplitude A [1]');
legend({'RF imageline','Low pass','High pass'});
ylim([	min([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25...
		max([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25]);

% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

%% ----------------------------
% 1D H-scan conv and rgb encoding
%-----------------------------
% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(imLineRF)));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
Hlo = conv(imLineRF, GHlo, 'same'); % Note: 'same' argument is important
Hhi = conv(imLineRF, GHhi, 'same');

% Envelope conversion (I was guessing that this is needed, please double-check)
Hlo = sqrt(abs(hilbert(Hlo)));
Hhi = sqrt(abs(hilbert(Hhi)));

Red_ImLine = Hlo./Hhi; % Stämmer det? Jag tror dessa varianter användes i någon av artiklarna
Blue_ImLine = Hhi./Hlo; % Stämmer det?

figure(3); clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channels for frame \newline',num2str(frame),' and line ',num2str(line)]);
legend({'Red channel','Blue channel'});

%% ---------------------------
% 2D H-scan conv and rgb encoding
%-----------------------------
clear RFlo RFhi
clear Bmodeslo Bmodeshi

% convolution
for k=1:128
	RFlo(:,k)=conv(RF(:,k),GHlo,'same')./sqrt(Elo);
	RFhi(:,k)=conv(RF(:,k),GHhi,'same')./sqrt(Ehi);
end

% envelope detection
Bmodeslo=sqrt(abs(hilbert(RFlo)));
Bmodeshi=sqrt(abs(hilbert(RFhi)));

%% plot H-scan filtering results
draw_pic(Bmodes(:,:), Bmodeslo, Bmodeshi, [] , 0.05,4);

%% rb colorcoding
BmodesrgbHr = zeros(shape(1),shape(2),3);
Bmodesrgb = zeros(shape(1),shape(2),3);
Bmodesrgblo = zeros(shape(1),shape(2),3);
Bmodesrgbhi = zeros(shape(1),shape(2),3);

BmodesrgbHr(:,:,1) = Bmodeshi(:,:)./mean2(Bmodeshi(:,:));
BmodesrgbHr(:,:,3) = Bmodeslo(:,:)./mean2(Bmodeslo(:,:));
Bmodesrgb(:,:,2) = Bmodes(:,:)./mean2(Bmodes(:,:));

Bmodesrgbhi(:,:,1,:) = BmodesrgbHr(:,:,1,:);
Bmodesrgblo(:,:,3,:) = BmodesrgbHr(:,:,3,:);

BmodesrgbHr = resRed(BmodesrgbHr);
Bmodesrgb = resRed(Bmodesrgb);
Bmodesrgbhi = resRed(Bmodesrgbhi);
Bmodesrgblo = resRed(Bmodesrgblo);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbHr,0.05,5);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes,BmodesrgbHr,0.05,6);

%% --------
% AT WORK
% --------
frame = 10; % Chose a frame, for dataset 1 one is sufficient for the report.
line=64; % Image line = 64 out or 128
load('C:\Users\Rebecca Viklund\Desktop\AMI project\AMI\Dataset_1\CG_contraction_1_omg2.mat');

% ----------------------
% Pre-processing of data
%-----------------------
RF = double(RF(1:2177,:,frame)); %1:2177 corresp to 4cm depth in this dataset..

% TGC - Only needed for dataset 1. TGC = "time gain compensation"
% It amplifies the signal proportional to the US depth since it weakens as
% it interacts with tissue.
TGC_VECT = linspace(1,10,size(RF,1))';
TGC_MAT = repmat(TGC_VECT,[1 size(RF,2)]);

RF = RF.*TGC_MAT;

% Here a "beamforming" is conducted. This is only needed for dataset 1.
% This is due to the collection method "plane-wave imaging" used which
% requires synthetic focusing which is provided by beamforming techniques.

BF = beamformDAQ(RF,128);
BF = BF - repmat(mean(BF,1),[size(BF,1) 1]);

% Use the beamformed version of RF for all computations below
RF = BF;

RF= single(squeeze(RF(:,:,:)));
Bmodes = single(sqrt(abs(hilbert(squeeze(RF(:,:,:))))));

% ----------------------
% H-scan 1b
%-----------------------
Fs = 35; % MHz Samplingfreq of the US not the signal
T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);

% GH low pass
b1 = 0.13;
ordlo = 8;
Hlo = hermiteH(ordlo, t./b1);
GHlo = exp(-(t./(b1)).^2).*Hlo;
GHlo = GHlo./sum(GHlo(:));
[pxxlo,flo] = pwelch(GHlo, hamming(512));
f_VECTlo = linspace(0,Fs/2,length(pxxlo));
p_NORMlo = 0.5*pxxlo./max(pxxlo);

% GH high pass
b2 = 0.135;
ordhi = 32;
Hhi = hermiteH(32, t./b2); % order 32
GHhi = exp(-(t./(b2)).^2).*Hhi;
GHhi = GHhi./sum(GHhi(:));
[pxxhi,fhi] = pwelch(GHhi,hamming(512));
f_VECThi = linspace(0,Fs/2,length(pxxhi));
p_NORMhi = 0.5*pxxhi./max(pxxhi);

% Plot one spectrum of one imageline for comparison with the GH spectra
imLineRF = double(squeeze(RF(:,line))); 
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));
p_NORM = sqrt(pxx)./max(sqrt(pxx));

figure(1);
subplot(2,1,1)
plot(t,GHlo,'b', 'LineWidth', 1.5);
title(['GH_{',num2str(ordlo),'}'])
ylim([min(GHlo,[],'all')*1.25 max(GHlo,[],'all')*1.25]);
subplot(2,1,2)
plot(t,GHhi,'r', 'LineWidth', 1.5);
title(['GH_{',num2str(ordhi),'}'])
ylim([min(GHhi,[],'all')*1.25 max(GHhi,[],'all')*1.25]);
sgtitle('Gaussian weighted Hermite polynomials')

figure(2); clf; hold on;
plot(f_VECT, p_NORM,'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECTlo, p_NORMlo, 'b', 'LineWidth', 1.5);
plot(f_VECThi, p_NORMhi, 'r', 'LineWidth', 1.5);
sgtitle(['PSD for frame ',num2str(frame),' and line ',num2str(line)]);
xlabel('Frequency f [MHz]'); ylabel('Amplitude A [1]');
legend({'RF imageline','Low pass','High pass'});
ylim([	min([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25...
		max([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25]);

% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

%% ----------------------------
% 1D H-scan conv and rgb encoding
%-----------------------------
% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(imLineRF)));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
Hlo = conv(imLineRF, GHlo, 'same'); % Note: 'same' argument is important
Hhi = conv(imLineRF, GHhi, 'same');

% Envelope conversion (I was guessing that this is needed, please double-check)
Hlo = sqrt(abs(hilbert(Hlo)));
Hhi = sqrt(abs(hilbert(Hhi)));

Red_ImLine = Hlo./Hhi; % Stämmer det? Jag tror dessa varianter användes i någon av artiklarna
Blue_ImLine = Hhi./Hlo; % Stämmer det?

figure(3); clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channels for frame \newline',num2str(frame),' and line ',num2str(line)]);
legend({'Red channel','Blue channel'});

%% ---------------------------
% 2D H-scan conv and rgb encoding
%-----------------------------
clear RFlo RFhi
clear Bmodeslo Bmodeshi

% convolution
for k=1:128
	RFlo(:,k)=conv(RF(:,k),GHlo,'same')./sqrt(Elo);
	RFhi(:,k)=conv(RF(:,k),GHhi,'same')./sqrt(Ehi);
end

% envelope detection
Bmodeslo=sqrt(abs(hilbert(RFlo)));
Bmodeshi=sqrt(abs(hilbert(RFhi)));

%% plot H-scan filtering results
draw_pic(Bmodes(:,:), Bmodeslo, Bmodeshi, [] , 0.05,4);

%% rb colorcoding
BmodesrgbH = zeros(shape(1),shape(2),3);
Bmodesrgb = zeros(shape(1),shape(2),3);
Bmodesrgblo = zeros(shape(1),shape(2),3);
Bmodesrgbhi = zeros(shape(1),shape(2),3);

BmodesrgbH(:,:,1) = Bmodeshi(:,:)./mean2(Bmodeshi(:,:));
BmodesrgbH(:,:,3) = Bmodeslo(:,:)./mean2(Bmodeslo(:,:));
Bmodesrgb(:,:,2) = Bmodes(:,:)./mean2(Bmodes(:,:));

Bmodesrgbhi(:,:,1,:) = BmodesrgbH(:,:,1,:);
Bmodesrgblo(:,:,3,:) = BmodesrgbH(:,:,3,:);

BmodesrgbH = resRed(BmodesrgbH);
Bmodesrgb = resRed(Bmodesrgb);
Bmodesrgbhi = resRed(Bmodesrgbhi);
Bmodesrgblo = resRed(Bmodesrgblo);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbH,0.05,5);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes,BmodesrgbH,0.05,6);

%% ---------------
% 1c Comparing muscle at rest and at work H-scan images
%-----------------
draw_pic2(BmodesrgbHr,BmodesrgbH,0.05,11);

%% load dataset 2 stimulated contraction w ~1.6Hz 1 muscle complex
% should be conducted on several frames and is to be filtered using
% svd only.
% ----------------------
% Pre-processing data
%-----------------------
rfmat_dsf = single(load('181023_1311_rs.mat').rfmat_downsampled);

Bmodes = sqrt(abs(hilbert(squeeze(rfmat_dsf(:,:,:)))));
shape = size(Bmodes);
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));

% ---------------------
% SVD 2 - Viktors part
%-----------------------
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(25:end, 25:end) = 0;
Snew(1:20,1:20) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

% MY COMPUTER CAN'T HANDLE THIS DATASET!

%% load dataset 3 voluntary contraction w ~5-25Hz 2-4 muscle complexes
% should be conducted on several frames and is to be filtered using both
% svd and H-scan
% ----------------------
% Pre-processing data
%-----------------------
load("RF_MAT.mat");
RF_MAT= single(squeeze(RF_MAT(:,:,:)));
Bmodes = single(sqrt(abs(hilbert(squeeze(RF_MAT(:,:,:))))));
shape = size(Bmodes);

% ---------------------
% SVD 3 - Viktors part
%-----------------------


% ---------------------
% H-scan 3
%-----------------------
% checking filter settings in 1D by taking a line/col of the img for
% freq analysis
Fs = 35; % MHz Samplingfreq of the US not the signal
T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);

% GH low pass
b1 = 0.13;
ordlo = 8;
Hlo = hermiteH(ordlo, t./b1);
GHlo = exp(-(t./(b1)).^2).*Hlo;
GHlo = GHlo./sum(GHlo(:));
[pxxlo,flo] = pwelch(GHlo, hamming(512));
f_VECTlo = linspace(0,Fs/2,length(pxxlo));
p_NORMlo = 0.5*pxxlo./max(pxxlo);

% GH high pass
b2 = 0.135;
ordhi = 32;
Hhi = hermiteH(32, t./b2); % order 32
GHhi = exp(-(t./(b2)).^2).*Hhi;
GHhi = GHhi./sum(GHhi(:));
[pxxhi,fhi] = pwelch(GHhi,hamming(512));
f_VECThi = linspace(0,Fs/2,length(pxxhi));
p_NORMhi = 0.5*pxxhi./max(pxxhi);

% Plot one spectrum of one imageline for comparison with the GH spectra
imLineRF = double(squeeze(RF_MAT(:,line,frame))); 
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));
p_NORM = sqrt(pxx)./max(sqrt(pxx));

figure(1);
subplot(2,1,1)
plot(t,GHlo,'b', 'LineWidth', 1.5);
title(['GH_{',num2str(ordlo),'}'])
ylim([min(GHlo,[],'all')*1.25 max(GHlo,[],'all')*1.25]);
subplot(2,1,2)
plot(t,GHhi,'r', 'LineWidth', 1.5);
title(['GH_{',num2str(ordhi),'}'])
ylim([min(GHhi,[],'all')*1.25 max(GHhi,[],'all')*1.25]);
sgtitle('Gaussian weighted Hermite polynomials')

figure(2); clf; hold on;
plot(f_VECT, p_NORM,'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECTlo, p_NORMlo, 'b', 'LineWidth', 1.5);
plot(f_VECThi, p_NORMhi, 'r', 'LineWidth', 1.5);
sgtitle(['PSD for frame ',num2str(frame),' and line ',num2str(line)]);
xlabel('Frequency f [MHz]'); ylabel('Amplitude A [1]');
legend({'RF imageline','Low pass','High pass'});
ylim([	min([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25...
		max([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25]);

% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

%% ----------------------------
% 1D H-scan conv and rgb encoding
%-----------------------------
% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(imLineRF)));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
Hlo = conv(imLineRF, GHlo, 'same'); % Note: 'same' argument is important
Hhi = conv(imLineRF, GHhi, 'same');

% Envelope conversion (I was guessing that this is needed, please double-check)
Hlo = sqrt(abs(hilbert(Hlo)));
Hhi = sqrt(abs(hilbert(Hhi)));

Red_ImLine = Hlo./Hhi; % Stämmer det? Jag tror dessa varianter användes i någon av artiklarna
Blue_ImLine = Hhi./Hlo; % Stämmer det?

figure(3); clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channels for frame \newline',num2str(frame),' and line ',num2str(line)]);
legend({'Red channel','Blue channel'});

%% ---------------------------
% 2D H-scan conv and rgb encoding
%-----------------------------
clear RF_MATlo RF_MAThi
clear Bmodeslo Bmodeshi
noframes=30;

% convolution
for j=1:noframes
	for k=1:128
		RF_MATlo(:,k,j)=conv(RF_MAT(:,k,j),GHlo,'same')./sqrt(Elo);
		RF_MAThi(:,k,j)=conv(RF_MAT(:,k,j),GHhi,'same')./sqrt(Ehi);
	end
end

% envelope detection
Bmodeslo=sqrt(abs(hilbert(RF_MATlo)));
Bmodeshi=sqrt(abs(hilbert(RF_MAThi)));

%% plot H-scan filtering results
draw_pic(Bmodes(:,:,1:noframes), Bmodeslo, Bmodeshi, [] , 0.05,9);

%% rb colorcoding
BmodesrgbH = zeros(shape(1),shape(2),3,noframes);
Bmodesrgb = zeros(shape(1),shape(2),3,noframes);
Bmodesrgblo = zeros(shape(1),shape(2),3,noframes);
Bmodesrgbhi = zeros(shape(1),shape(2),3,noframes);

for j=1:noframes
	BmodesrgbH(:,:,1,j) = squeeze(Bmodeshi(:,:,j))./mean2(squeeze(Bmodeshi(:,:,j)));
	BmodesrgbH(:,:,3,j) = squeeze(Bmodeslo(:,:,j))./mean2(squeeze(Bmodeslo(:,:,j)));
 	Bmodesrgb(:,:,2,j) = squeeze(Bmodes(:,:,j))./mean2(squeeze(Bmodes(:,:,j)));
end

Bmodesrgbhi(:,:,1,:)=BmodesrgbH(:,:,1,:);
Bmodesrgblo(:,:,3,:)=BmodesrgbH(:,:,3,:);

for j=1:noframes
    BmodesrgbH(:,:,:,j) = resRed(BmodesrgbH(:,:,:,j));
    Bmodesrgb(:,:,:,j) = resRed(Bmodesrgb(:,:,:,j));
    Bmodesrgbhi(:,:,:,j) = resRed(Bmodesrgbhi(:,:,:,j));
    Bmodesrgblo(:,:,:,j) = resRed(Bmodesrgblo(:,:,:,j));
end

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbH,0.05,10);

%%
% plot comparison between B-mode and H-scan
draw_pic2(Bmodes(:,:,1:noframes),BmodesrgbH,0.1,11);

%%
% plot a frame between 1 and noframes
frame1=18;
draw_pic2(Bmodes(:,:,frame1),BmodesrgbH(:,:,:,frame1),0.1,11);
%% load TVI (Validation) data
% ----------------------
% Pre-processing data
%-----------------------
fn_STR = '181023_1311tvi_rs.mat';
h = matfile(fn_STR);
tvif_d = single(h.tvi_downsampled(:,:,1:250));
tvif_d = sqrt(abs(hilbert(squeeze(tvif_d(:,:,:)))));

%% SVD vs TVI 
% Prefiltering
% median3 and butter order 4 on both TVI and SVD, keeping the freqs in
% range 5-25 Hz

% tripple median filtering
Bmodes_fin=medfilt3(Bmodes_new);
tvi_fin=medfilt3(tvif);

% spatial filtering
C = conv2(A, B8);

% freq filtering
% idea: bandpass on every pixel intensity over activation time 5-15(30) Hz
% for dataset 2 but around 1.5 Hz for dataset 3
[b,a] = butter(4,[5 25]./(Fsamp/2),'bandpass');
y_filt = filtfilt(b,a,y);

%% H-scan vs TVI
% Prefiltering
% median3 and butter order 4 on both TVI and SVD, keeping the freqs in
% range 5-25 Hz

% tripple median filtering
Bmodes_fin=medfilt3(Bmodes_new);
tvi_fin=medfilt3(tvif);

% spatial filtering
C = conv2(A, B8);

% freq filtering
% idea: bandpass on every pixelintensity over activationtime 5-15(30) Hz
% for dataset 2 but around 1.5 Hz for dataset 3
[b1,b2] = butter(4,[5 25]./(Fsamp/2),'bandpass');
y_filt = filtfilt(b1,b2,y);

%%
function filteredI = resRed(I)
% Reduces resolution on RGB channels to improve visibility of changes in
% H-scan imaging by taking the median of every 3x3 area and overwriting all
% pixels within it to its value.
% Input: RGB image I 
% Output: Down resoulutioned version of the RGB image.
if size(I,4)==1
    for j=2:3:(size(I,1)-1)
        for k=2:3:(size(I,2)-1)
            Imed=median([I(j-1,k-1:k+1,:,:) ...
                I(j,k-1:k+1,:,:) I(j+1,k-1:k+1,:,:)]);
                ImedMat=repmat(Imed,[3,3,1,1]);
                I(j-1:j+1,k-1:k+1,:,:)=ImedMat;
        end
    end
    filteredI=I;
else
    for j=2:3:(size(I,1)-1)
        for k=2:3:(size(I,2)-1)
            Imed=median([I(j-1,k-1:k+1,:) ...
                I(j,k-1:k+1,:) I(j+1,k-1:k+1,:)]);
                ImedMat=repmat(Imed,[3,3,1]);
                I(j-1:j+1,k-1:k+1,:)=ImedMat;
        end
    end
    filteredI=I;
end
end

function [filteredI,mask] = maskFilter(I,thresh)
	% Get a mask using locally adaptive thresholding.
	mask = imbinarize(I, thresh);
	% Dilate it a bit to make the spots bigger.
	mask = imdilate(mask,strel("disk",r,n));
	% Repair the image by using regionfill() 
	% to smear in the surrounding pixels values into the white spots.
	repairedImage = regionfill(I, mask);
	filteredI=repairedImage;
end

function draw_pic(mat1, mat2, mat3, mat4, delay, fignum)
% plots 3-4 Bmode images with specified delay and option of including
% figure number
	lowi = 0.0;
	highi = 0.1;

	if fignum
		figure(fignum);
	else 
		figure
	end
	
	% rgb or not condition
	if size(mat2,3)==3
		Q = size(mat2,4);
		
		W1 = mat1(:,:,:,1);
		%W1 = mat2gray(W1);

    	W2 = mat2(:,:,:,1);

		W3 = mat3(:,:,:,1);

		W4 = mat4(:,:,:,1);

		% Maximize window.
		g = gcf;
% 		g.WindowState = 'maximized';
		drawnow;
	
    	subplot(2,2,1);
    	img1 = imagesc(W1, [lowi, highi]);
		title('Original Bmode color channel')
	
    	subplot(2,2,2);
		img2 = imagesc(W2, [lowi, highi]);
		title('Low pass filtered Bmodes color channel')

		subplot(2,2,3);
		img3 = imagesc(W3, [lowi, highi]);
		title('High pass filtered Bmodes color channel')

		subplot(2,2,4);
		img4 = imagesc(W4, [lowi, highi]);
		title('Rgb coded bandpass filtered Bmode')
		sgtitle('H-scan images compared to the original B-modes');

		pause(2);

		for K = 2:Q
        W1 = mat1(:,:,:,K);
		%W1 = mat2gray(W1);
        W1 = imadjust(W1, [0,1],[]);
		
        W2 = mat2(:,:,:,K);
		W2 = imadjust(W2, [0,1],[]);

		W3 = mat3(:,:,:,K);
        W3 = imadjust(W3, [0,1],[]);

		W4 = mat4(:,:,:,K);
		W4 = imadjust(W4, [0,1],[]);

        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
		set(img3, 'CData', W3);
		set(img4, 'CData', W4);
		
        %drawnow limitrate;
        drawnow();
        clim([0,1])
        pause(delay);
        disp(K);
		end
		
	else
		Q = size(mat1,3);
		
		W1 = mat1(:,:,1);
		W1 = mat2gray(W1);
	
    	W2 = mat2(:,:,1);
		W2 = mat2gray(W2);

		W3 = mat3(:,:,1);
		W3 = mat2gray(W3);

		% Maximize window.
		g = gcf;
    	g.WindowState = 'maximized';
		drawnow;
	
    	subplot(1,3,1);
    	img1 = imagesc(W1, [lowi, highi]);colormap gray;
		title('Original Bmode')
	
    	subplot(1,3,2);
		img2 = imagesc(W2, [lowi, highi]);
		title('Low pass filtered Bmode')

		subplot(1,3,3);
		img3 = imagesc(W3, [lowi, highi]);
		title('High pass filtered Bmode')

		sgtitle('H-scan images compared to the original B-modes');

		pause(2);
		for K = 2:Q
        W1 = mat1(:,:,K);
		W1 = mat2gray(W1);
        W1 = imadjust(W1, [0,1],[]);

        W2 = mat2(:,:,K);
        W2 = mat2gray(W2);
		W2 = imadjust(W2, [0,1],[]);

		W3 = mat3(:,:,K);
		W3 = mat2gray(W3);
        W3 = imadjust(W3, [0,1],[]);

        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
		set(img3, 'CData', W3);

        %drawnow limitrate;
        drawnow();
        clim([0,1])
        pause(delay);
        disp(K);
		end
	end
end

function draw_pic2(mat1, mat2, delay, fignum)
% plots 3-4 Bmode images with specified delay and option of including
% figure number
	lowi = 0.0;
	highi = 0.1;

	if fignum
		figure(fignum);
	else 
		figure
	end
		Q = size(mat1,3);
		
	if size(mat1,3)==3
		W1 = mat1(:,:,:,1);
	else
		W1 = mat1(:,:,1);
		W1 = mat2gray(W1);
	end
    	W2 = mat2(:,:,:,1);

		% Maximize window.
		g = gcf;
% 		g.WindowState = 'maximized';
		drawnow;
	
    	subplot(1,2,1);
    	img1 = imagesc(W1, [lowi, highi]);colormap gray;
        if size(mat1,3)==3
		    title('At rest')
        else
            title('B-mode')
        end
    	subplot(1,2,2);
		img2 = imagesc(W2, [lowi, highi]);
		if size(mat1,3)==3
        	title('Contracted')
            sgtitle(['Comparison of activity between fully contracted and ' ...
                'at rest'])
        else
            title('H-scan')
        end
		pause(2);

		for K = 2:Q
		if size(mat1,3)==3
		W1 = mat1(:,:,:,K);
		else
        W1 = mat1(:,:,K);
		W1 = mat2gray(W1);
		end
        W1 = imadjust(W1, [0,1],[]);
		
        W2 = mat2(:,:,:,K);
		W2 = imadjust(W2, [0,1],[]);

        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
		
        %drawnow limitrate;
        drawnow();
        clim([0,.95])
        pause(delay);
        disp(K);
		end
end