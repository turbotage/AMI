%  AMI Project 2022

%% load dataset 1 ... only image contraction w ~ Hz full activation vs full rest
% should be conducted on one frame with only the H-scan for both full
% activation and at rest.

frame = 10; % Chose a frame, for dataset 1 one is sufficient for the report.
line=64; % Image line = 64 out or 128

%% at rest
load('C:\Users\revi0014\Desktop\AMI\Dataset_1\CG_rest_2_omg2.mat');

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

% ----------------------
% H-scan 1a
%-----------------------
Fs = 35; %MHz Samplingfreq of the US not the signal
T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);

b1 = 0.08;
H2 = (4.*(t./b1).^2 - 2);
GH2 = exp(-(t./(b1)).^2).*H2;
GH2 = GH2./sum(GH2(:));
[pxx2,f2] = pwelch(GH2, hamming(512));
f_VECT2 = linspace(0,Fs/2,length(pxx2));

b2 = 0.07;
H8 = (256.*(t./b2).^8 - 3584.*(t./b2).^6 + 13440.*(t./b2).^4 - 13440.*(t./b2).^2 + 1680);
GH8 = exp(-(t./(b2)).^2).*H8;
GH8 = GH8./sum(GH8(:));
[pxx8,f1] = pwelch(GH8,hamming(512));
f_VECT8 = linspace(0,Fs/2,length(pxx8));

% Plot one spectrum of one imageline for comparison with the GH spectra
imLineRF = double(squeeze(RF(:,line))); 
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));

figure(1);
subplot(2,1,1); plot(t, GH2, 'b'); ylabel('GH_2(t)');
subplot(2,1,2); plot(t, GH8, 'r'); ylabel('GH_8(t)');xlabel('Time, us');

figure(2); clf; hold on;
plot(f_VECT, sqrt(pxx)./max(sqrt(pxx)),'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECT2, 0.5*pxx2./max(pxx2), 'b', 'LineWidth', 1.5);
plot(f_VECT8, 0.5*pxx8./max(pxx8), 'r', 'LineWidth', 1.5);
xlabel('Frequency, MHz');
ylabel('Amplitude, n.u.');
legend({'RF imageline','GH_2','GH_8'});
hold off;

% ----------------------------
% H-scan conv and rgb encoding
%-----------------------------
% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(imLineRF)));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
H2 = conv(imLineRF, GH2, 'same'); % Note: 'same' argument is important
H8 = conv(imLineRF, GH8, 'same');

% Envelope conversion (I was guessing that this is needed, please double-check)
H2 = sqrt(abs(hilbert(H2)));
H8 = sqrt(abs(hilbert(H8)));

Red_ImLine = H2./H8; % Stämmer det? Jag tror dessa varianter användes i någon av artiklarna
Blue_ImLine = H8./H2; % Stämmer det?

figure(3); clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
legend({'Red channel','Blue channel'});

%% during contraction
load('C:\Users\revi0014\Desktop\AMI\Dataset_1\CG_contraction_1_omg2.mat');

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

% ----------------------
% H-scan 1b
%-----------------------
Fs = 35; %MHz Samplingfreq of the US not the signal
T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);

b1 = 0.08;
H2 = (4.*(t./b1).^2 - 2);
GH2 = exp(-(t./(b1)).^2).*H2;
GH2 = GH2./sum(GH2(:));
[pxx2,f2] = pwelch(GH2, hamming(512));
f_VECT2 = linspace(0,Fs/2,length(pxx2));

b2 = 0.07;
H8 = (256.*(t./b2).^8 - 3584.*(t./b2).^6 + 13440.*(t./b2).^4 - 13440.*(t./b2).^2 + 1680);
GH8 = exp(-(t./(b2)).^2).*H8;
GH8 = GH8./sum(GH8(:));
[pxx8,f1] = pwelch(GH8,hamming(512));
f_VECT8 = linspace(0,Fs/2,length(pxx8));

% Plot one spectrum of one imageline for comparison with the GH spectra
imLineRF = double(squeeze(RF(:,line))); 
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));

figure(4);
subplot(2,1,1); plot(t, GH2, 'b'); ylabel('GH_2(t)');
subplot(2,1,2); plot(t, GH8, 'r'); ylabel('GH_8(t)');xlabel('Time, us');

figure(5); clf; hold on;
plot(f_VECT, sqrt(pxx)./max(sqrt(pxx)),'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECT2, 0.5*pxx2./max(pxx2), 'b', 'LineWidth', 1.5);
plot(f_VECT8, 0.5*pxx8./max(pxx8), 'r', 'LineWidth', 1.5);
xlabel('Frequency, MHz');
ylabel('Amplitude, n.u.');
legend({'RF imageline','GH_2','GH_8'});
hold off;

% ----------------------------
% H-scan conv and rgb encoding
%-----------------------------
% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(imLineRF)));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
H2 = conv(imLineRF, GH2, 'same'); % Note: 'same' argument is important
H8 = conv(imLineRF, GH8, 'same');

% Envelope conversion (I was guessing that this is needed, please double-check)
H2 = sqrt(abs(hilbert(H2)));
H8 = sqrt(abs(hilbert(H8)));

Red_ImLine = H2./H8; % Stämmer det? Jag tror dessa varianter användes i någon av artiklarna
Blue_ImLine = H8./H2; % Stämmer det?

figure(6); clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
legend({'Red channel','Blue channel'});

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
clear Bmodes

%% ---------------------
% SVD 2
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
RF_MAT_f=reshape(RF_MAT, shape(1)*shape(2), shape(3));

%% ---------------------
% SVD 3
%-----------------------
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(25:end, 25:end) = 0;
Snew(1:20,1:20) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%% ---------------------
% H-scan 3
%-----------------------
b1 = .1;
b2 = .07;
fs= 35; %should be for the ultrasound, not 1000;
T_duration = 50; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*fs);

GH2=exp(-(t/b1).^2).*(4*(t./b1).^2-2);
GH2 = GH2./sum(GH2(:));
[pxx2,f2] = pwelch(GH2, hamming(512));
f_VECT2 = linspace(0,fs/2,length(pxx2));
p_NORM2 = 0.5*pxx2./max(pxx2);

GH8=exp(-(t/b2).^2).*((2*(t./b2)).^8-3584*(t./b2).^6+13440*((t./b2).^4-(t./b2).^2)+1680);
GH8 = GH8./sum(GH8(:));
[pxx8,f1] = pwelch(GH8,hamming(512));
f_VECT8 = linspace(0,fs/2,length(pxx8));
p_NORM8 = 0.5*pxx8./max(pxx8);

E2 = 1*3*sqrt(pi/2);
E8 = 1*3*5*7*9*11*13*15*sqrt(pi/2);

sgtitle('Gaussian weighted Hermite polynomials')
subplot(2,1,1)
plot(t,GH2,'b', 'LineWidth', 1.5);
title('GH_2')
ylim([min(GH2,[],'all')*1.25 max(GH2,[],'all')*1.25]);
subplot(2,1,2)
plot(t,GH8,'r', 'LineWidth', 1.5);
title('GH_8')
ylim([min(GH8,[],'all')*1.25 max(GH8,[],'all')*1.25]);

% Plot one spectrum of one imageline for comparison with the GH spectra
frame=45;	% out of 500
line=64;	% out of 128

imLineRF = double(squeeze(RF_MAT(:,line,frame)));
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,fs/2,length(f));
p_NORM = sqrt(pxx)./max(sqrt(pxx));

figure; clf; hold on;
plot(f_VECT, p_NORM,'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECT2, p_NORM2, 'b', 'LineWidth', 1.5);
plot(f_VECT8, p_NORM8, 'r', 'LineWidth', 1.5);
sgtitle(['PSD for frame ',num2str(frame),' and line ',num2str(line)]);
xlabel('Frequency f [MHz]'); ylabel('Amplitude A [1]');
legend({'RF imageline','Low pass (GH2)','High pass (GH8)'});
ylim([	min([p_NORM2 p_NORM8 p_NORM],[],'all')*1.25...
		max([p_NORM2 p_NORM8 p_NORM],[],'all')*1.25]);

% ---------------------------
% H-scan conv and rgb encoding
%-----------------------------
% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(RF_MAT(:,line,frame))));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
H2 = conv(RF_MAT(:,line,frame), GH2, 'same');
H8 = conv(RF_MAT(:,line,frame), GH8, 'same');
% envelope conversion (I was guessing that this is needed, please double-check)

H2 = sqrt(abs(hilbert(H2)));
H8 = sqrt(abs(hilbert(H8)));

Red_ImLine = H8./H2; % Stämmer det? Jag tror dessa varianter användes i någon av artiklarna
Blue_ImLine = H2./H8; % Stämmer det?

figure; clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channels for frame \newline',num2str(frame),' and line ',num2str(line)]);
legend({'Red channel','Blue channel'});

%% ---------------------------
% 2D H-scan conv and rgb encoding
%-----------------------------
clear RF_MAT2 RF_MAT8
clear Bmodes2 Bmodes8

% convolution
noframes=50;
for j=1:noframes
	for k=1:128
		RF_MAT2(:,k,j)=conv(RF_MAT(:,k,j),GH2,'same')./sqrt(E2);
		RF_MAT8(:,k,j)=conv(RF_MAT(:,k,j),GH8,'same')./sqrt(E8);
	end
end

% envelope detection
Bmodes2=double(sqrt(abs(hilbert(RF_MAT2))));
Bmodes8=sqrt(abs(hilbert(RF_MAT8)));

% plot H-scan filtering results
draw_pic(Bmodes(:,:,1:noframes), Bmodes2, Bmodes8, [] , 0.05);

%% Colorcode in 2D
% Combination with colorcoding in freq space to RGB and then combined into 
% a new image RGB colormap R for low freqs, B for high, G for original sig
% envelope.

% RGB channeling intensites G=signal, B=GH8, R=GH2
% mymap = [	1 0 0 ;
%			0 1 0 ;
%			0 0 1];


Bmodesrgb = zeros(shape(1),shape(2),3,noframes);
Bmodes2rgb = zeros(shape(1),shape(2),3,noframes);
Bmodes8rgb = zeros(shape(1),shape(2),3,noframes);
Bmodestotrgb = zeros(shape(1),shape(2),3,noframes);

for m=1:noframes
	Bmodesrgb(:,:,:,m)=ind2rgb(squeeze(Bmodes(:,:,m)),'turbo');
	Bmodes2rgb(:,:,:,m)=ind2rgb(squeeze(Bmodes2(:,:,m)),'turbo');
	Bmodes8rgb(:,:,:,m)=ind2rgb(squeeze(Bmodes8(:,:,m)),'turbo');
end

Bmodestotrgb(:,:,2,:)=Bmodesrgb(:,:,2,:);
Bmodestotrgb(:,:,1,:)=Bmodes8rgb(:,:,1,:);
Bmodestotrgb(:,:,3,:)=Bmodes2rgb(:,:,3,:);

draw_pic2(Bmodesrgb(:,:,:,:), Bmodestotrgb(:,:,:,:));

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
function [filteredI,mask]=maskFilter(I,thresh)
	% Get a mask using locally adaptive thresholding.
	mask = imbinarize(I, thresh);
	% Dilate it a bit to make the spots bigger.
	mask = imdilate(mask,strel("disk",r,n));
	% Repair the image by using regionfill() 
	% to smear in the surrounding pixels values into the white spots.
	repairedImage = regionfill(I, mask);
	filteredI=repairedImage;
end

function draw_pic(mat1, mat2, mat3, mat4, delay)
    lowi = 0.0;
    highi = 0.1;

    figure(1);
    Q = size(mat1,3);

    W1 = mat1(:,:,1);
	W1 = mat2gray(W1);

    W2 = mat2(:,:,1);
	W2 = mat2gray(W2);

	W3 = mat3(:,:,1);
	W3 = mat2gray(W3);

	if mat4
		W4 = mat4(:,:,1);
		W4 = mat2gray(W4);
	end

	% Maximize window.
	g = gcf;
	g.WindowState = 'maximized';
	drawnow;

    subplot(2,2,1);
	W1 = mat2gray(W1);
    img1 = imagesc(W1, [lowi, highi]),colormap gray;

    subplot(2,2,2);
	W2 = mat2gray(W2);
    img2 = imagesc(W2, [lowi, highi]);

	subplot(2,2,3);
	W3 = mat2gray(W3);
    img3 = imagesc(W3, [lowi, highi]);
	

	if mat4
		subplot(2,2,4);
		W4 = mat2gray(W4);
		img4 = imagesc(W4, [lowi, highi]);
	end

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

		if mat4
			W4 = mat4(:,:,K);
			W4 = mat2gray(W4);
			W4 = imadjust(W4, [0,1],[]);
		end

        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
		set(img3, 'CData', W3);
		if mat4
			set(img3, 'CData', W4);
		end
        %drawnow limitrate;
        drawnow();
        caxis([0,1])
        pause(delay);
        disp(K);
    end
end

function draw_pic2(mat1, mat2)
    lowi = 0.0;
    highi = 0.1;

    figure(1);
    Q = size(mat1,4);
    W1 = mat1(:,:,:,1);
    %W1 = imadjust(W1, [0,1],[]);

    W2 = mat2(:,:,:,2);
    %W2 = imadjust(W2, [0,1],[]);


    subplot(1,2,1);
    img1 = imagesc(W1, [lowi, highi]);

    subplot(1,2,2);
    img2 = imagesc(W2, [lowi, highi]);
	colormap turbo
	colorbar;
    
    pause(2);
    for K = 2:Q
        W1 = mat1(:,:,:,K);
        W1 = imadjust(W1, [0,1],[]);
        W2 = mat2(:,:,:,K);
        W2 = imadjust(W2, [0,1],[]);
        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
        %drawnow limitrate;
        drawnow();
        caxis([0,1])
        pause(0.02);
        disp(K);
    end
end
