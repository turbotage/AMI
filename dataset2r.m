% AMI Project 2022 - Rebecca
% Dataset 2 - H-scan
%--------------------------------------------------------------------------
% load dataset 2 stimulated contraction w ~1.6Hz i.e. 2-3 Hz on 1 muscle complex
% should be conducted on several frames and is to be filtered using
% svd and H-scan?
% ----------------------
% Pre-processing data
%-----------------------
% 1311 - 8Hz stim
rfmat = single(load('181023_1311_rs.mat').rfmat_downsampled);
rfmat = rfmat(1:1000,:,:);
TGC_VECT = linspace(1,10,size(rfmat,1))';
TCG_MAT = repmat(TGC_VECT,[1 size(rfmat,2)]);
rfmat=rfmat.*TCG_MAT;

roimat = load('181023_1311_ROI.mat');
roimat = roimat.ROI;

Bmodes = sqrt(abs(hilbert(squeeze(rfmat(:,:,:)))));
shape = size(Bmodes);

%% ---------------------
% SVD vs TVI - Viktors part
%-----------------------
% --------------------
% Pre-processing TVI data 
%---------------------
% fn_STR = '181023_1311tvi_rs.mat';
% h = matfile(fn_STR);
% tvif_d = single(h.tvi_downsampled(:,:,1:250));
% tvif_d = sqrt(abs(hilbert(squeeze(tvif_d(:,:,:)))));
%-------------------------
% Filtering and comparing
%-------------------------
% Prefiltering
% median3 and butter order 4 on both TVI and SVD, keeping the freqs in
% range 5-25 Hz

% tripple median filtering
% Bmodes_fin=medfilt3(Bmodes_new);
% tvi_fin=medfilt3(tvif);

% spatial filtering
% C = conv2(A, B8);

% freq filtering
% idea: bandpass on every pixel intensity over activation time 5-15(30) Hz
% for dataset 2 but around 1.5 Hz for dataset 3
% [b,a] = butter(4,[5 25]./(Fsamp/2),'bandpass');
% y_filt = filtfilt(b,a,y);

% ---------------------
% H-scan 3
%-----------------------
% checking filter settings in 1D by taking a line/col of the img for
% freq analysis
Fs = 40; % MHz Samplingfreq of the US not the signal
T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);

% GH low pass
b1 = 0.340;
ordlo = 2;
Hlo = hermiteH(ordlo, t./b1);
GHlo = exp(-(t./(b1)).^2).*Hlo;
GHlo = GHlo./sum(GHlo(:));

% GH high pass
b2 = 0.135;
ordhi = 32;
Hhi = hermiteH(32, t./b2); % order 32
GHhi = exp(-(t./(b2)).^2).*Hhi;
GHhi = GHhi./sum(GHhi(:));

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

%% Plot one spectrum of one imageline for comparison with the GH spectra
frame=floor(size(rfmat,3)/2); % 
line=floor(size(rfmat,2)/2); %

[pxxlo,flo] = pwelch(GHlo, hamming(512));
f_VECTlo = linspace(0,Fs/2,length(flo));
p_NORMlo = 0.5*pxxlo./max(pxxlo);

[pxxhi,fhi] = pwelch(GHhi,hamming(512));
f_VECThi = linspace(0,Fs/2,length(fhi));
p_NORMhi = 0.5*pxxhi./max(pxxhi);

imLineRF = double(squeeze(rfmat(:,line,frame))); 
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));
p_NORM = sqrt(pxx)./max(sqrt(pxx));

figure(2); clf; hold on;
plot(f_VECT, p_NORM,'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECTlo, p_NORMlo, 'b', 'LineWidth', 1.5);
plot(f_VECThi, p_NORMhi, 'r', 'LineWidth', 1.5);
sgtitle(['PSD for frame ',num2str(frame),' and line ',num2str(line)]);
xlabel('Frequency f [MHz]'); ylabel('Amplitude A [1]');
legend({'RF imageline','Low pass','High pass'});
ylim([	min([p_NORMlo p_NORMhi [p_NORM;ones(length(p_NORMlo)-length(p_NORM),1)]],[],'all')*1.15...
		max([p_NORMlo p_NORMhi [p_NORM;zeros(length(p_NORMlo)-length(p_NORM),1)]],[],'all')*1.15]);

% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

%% ---------------------------
% 2D H-scan conv and rgb encoding
%-----------------------------
clear RF_MATlo RF_MAThi
clear Bmodeslo Bmodeshi
noframes=size(rfmat,3);
frames=1:floor(noframes);

% convolution
for j=1:noframes
	for k=1:128
		RF_MATlo(:,k,j)=conv(rfmat(:,k,j),GHlo,'same')./sqrt(Elo);
		RF_MAThi(:,k,j)=conv(rfmat(:,k,j),GHhi,'same')./sqrt(Ehi);
	end
end

% envelope detection
Bmodeslo=sqrt(abs(hilbert(RF_MATlo)));
Bmodeshi=sqrt(abs(hilbert(RF_MAThi)));

% rb colorcoding
BmodesrgbH = zeros(shape(1),shape(2),3,noframes);
Bmodesrgb = zeros(shape(1),shape(2),3,noframes);
Bmodesrgblo = zeros(shape(1),shape(2),3,noframes);
Bmodesrgbhi = zeros(shape(1),shape(2),3,noframes);

for j=1:floor(noframes)
	BmodesrgbH(:,:,1,j) = (Bmodeshi(:,:,j)-min(Bmodeshi(:,:,j),[],'all')) / (max(Bmodeshi(:,:,j),[],'all')-min(Bmodeshi(:,:,j),[],'all')); %Bmodeshi(:,:)./mean2(Bmodeshi(:,:));
	BmodesrgbH(:,:,3,j) = (Bmodeslo(:,:,j)-min(Bmodeslo(:,:,j),[],'all')) / (max(Bmodeslo(:,:,j),[],'all')-min(Bmodeslo(:,:,j),[],'all'));
 	Bmodesrgb(:,:,2,j) = (Bmodes(:,:,j)-min(Bmodes(:,:,j),[],'all')) / (max(Bmodes(:,:,j),[],'all')-min(Bmodes(:,:,j),[],'all'));
end

[Bmodesrgbhi(:,:,1,:),~,Bmodesrgblo(:,:,3,:)]=imsplit(BmodesrgbH);

for j=1:floor(noframes)
    BmodesrgbH(:,:,1,j) = medfilt2(BmodesrgbH(:,:,1,j));
    BmodesrgbH(:,:,3,j) = medfilt2(BmodesrgbH(:,:,3,j));
    Bmodesrgb(:,:,2,j) = medfilt2(Bmodesrgb(:,:,2,j));
    Bmodesrgbhi(:,:,1,j) = medfilt2(Bmodesrgbhi(:,:,1,j));
    Bmodesrgblo(:,:,3,j) = medfilt2(Bmodesrgblo(:,:,3,j));
end

%% plot H-scan filtering results
draw_pic(Bmodes(:,:,1:noframes), Bmodeslo, Bmodeshi, [] , 0.05,9);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbH,0.05,10);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes(:,:,1:floor(noframes/4)),BmodesrgbH,0.1,11);

%% Mask the image
for j=1:noframes
    %Hmasked(:,:,:,j) = bsxfun(@times, BmodesrgbH(:,:,:,j), cast(roimat, 'like', BmodesrgbH(:,:,:,j)));
    Hmasked(:,:,:,j) = BmodesrgbH(:,:,:,j).*repmat(roimat,[1,1,3]);
end
%% plot mask vs H-scan
opacity=.1;
figure
W1=BmodesrgbH(:,:,:,1);
W2=roimat;
W3=Hmasked(:,:,:,1);

subplot(1,2,1)
img1=image(W1,'CDataMapping','scaled');
hold on;
img2=image(W2,'CDataMapping','scaled','AlphaData',opacity*ones(size(roimat)));colormap gray;
title(['H-scan with ROI overlapped at opacity ',num2str(opacity)])
subplot(1,2,2)
img3=image(W3,'CDataMapping','scaled');
title('masked out image area');
pause(2)
for j=1:floor(noframes/4)
    W1 = BmodesrgbH(:,:,:,j);
    W1 = imadjust(W1, [0,1],[]);

    W3 = Hmasked(:,:,:,j);
    W3 = imadjust(W3, [0,1],[]);
    
    set(img1, 'CData', W1);
    set(img2, 'CData', W2);
    set(img3, 'CData', W3);

    %drawnow limitrate;
    drawnow();
    clim([0,1])
    pause(.05);
    disp(j);
end

%% analyze intensity in a pixel; over measurement time (frames)
depth=floor(size(BmodesrgbH,1)/2);

for j=1:noframes
    pixelHred(j)=Bmodesrgbhi(depth,line,1,j);
    pixelHblu(j)=Bmodesrgblo(depth,line,3,j);
end

pixelHred(isnan(pixelHred))=0;
pixelHblu(isnan(pixelHblu))=0;

% pixelHred=pixelHred/pixelHblu;
% pixelHred=pixelHblu/pixelHred;

[pks,locs] = findpeaks(pixelHred,frames);
[pks1,locs1] = findpeaks(pixelHblu,frames);

pixelHredMax=makima(locs,pks,frames);
pixelHbluMax=makima(locs1,pks1,frames);

figure; clf; hold on;
plot(frames,pixelHred,'-r');
plot(frames,pixelHblu','-b');
plot(frames,pixelHredMax,'-r','LineWidth',1.5);
plot(frames,pixelHbluMax,'-b','LineWidth',1.5);
%ylim([0 1]);
xlabel('Frames t [ms]'); ylabel('Channel intensity I [%]');
sgtitle(['Comparison of channel intensities in pixel (',num2str(line),'x',num2str(depth),')']);
legend({'Red channel','Blue channel','Red trend','Blue trend'});

%% analyze std of intensity in each frame; over measurement time (frames)

for j=1:noframes
    frameHred(:,j) = std2(Bmodesrgbhi(:,:,1,j));
    frameHblu(:,j) = std2(Bmodesrgblo(:,:,3,j));
end

frameHred(isnan(frameHred))=0;
frameHblu(isnan(frameHblu))=0;

% frameHred=frameHred/frameHblu;
% frameHred=frameHblu/frameHred;

figure; clf; hold on;
plot(frames,frameHred,'-r');
plot(frames,frameHblu','-b');
% ylim([0 1]);
xlabel('Frames t [ms]'); ylabel('Channel intensity std \sigma [1]');
sgtitle(['Comparison of channel std of intensities in frames']);
legend({'Red channel','Blue channel','Red trend','Blue trend'});

%% analyze intensity in a line; over signal time (depth)
nodepths=size(BmodesrgbH,1);
depths=1:nodepths;

for j=1:nodepths
    lineHred(j)=Bmodesrgbhi(j,line,1,frame);
    lineHblu(j)=Bmodesrgblo(j,line,3,frame);
end

lineHred(isnan(lineHred))=0;
lineHblu(isnan(lineHblu))=0;

% lineHred=lineHred/lineHblu;
% lineHred=lineHblu/lineHred;

figure(3); clf; hold on;
plot(1:nodepths,lineHred,'r');
plot(1:nodepths,lineHblu,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channel intensities for frame ',num2str(frame),' over line ',num2str(line)]);
legend({'Red channel','Blue channel'});
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'Red channel','Blue channel','Red trend','Blue trend'});

%% analyze intensity in one frame; over signal time (depth)
nolines=size(BmodesrgbH,2);
lines=1:nolines;

for k=1:nodepths
    for j=1:nolines
        Tempr(j)=Bmodesrgbhi(k,j,1,frame);
        Tempb(j)=Bmodesrgblo(k,j,3,frame);
    end
    lineFrameHred(k)=mean(Tempr);
    lineFrameHblu(k)=mean(Tempb);
end

lineFrameHred(isnan(lineFrameHred))=0;
lineFrameHblu(isnan(lineFrameHblu))=0;

% lineFrameHred=lineFrameHred/lineFrameHblu;
% lineFrameHred=lineFrameHblu/lineFrameHred;

figure(3); clf; hold on;
plot(depths,lineFrameHred,'r');
plot(depths,lineFrameHblu,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channel mean intensities for\newlineframe ',num2str(frame),' over depth']);
legend({'Red channel','Blue channel'});
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'Red channel','Blue channel','Red trend','Blue trend'});

%% analyze mean intensity off all frames; over signal time (depth)

for l=1:noframes
    for k=1:nodepths
        for j=1:nolines
            Tempr(j)=Bmodesrgbhi(k,j,1,l);
            Tempb(j)=Bmodesrgblo(k,j,3,l);
        end
    lineFramesHred(k,l)=mean(Tempr);
    lineFramesHblu(k,l)=mean(Tempb);
    end
end

lineFramesHred(isnan(lineFramesHred))=0;
lineFramesHblu(isnan(lineFramesHblu))=0;

% lineFramesHred=lineFramesHred/lineFramesHblu;
% lineFramesHred=lineFramesHblu/lineFramesHred;

[X,Y] = meshgrid(frames,depths);
figure
g = gcf;
g.WindowState = 'maximized';
subplot(1,2,1)
surf(X,Y,lineFramesHred,'lineStyle','none');
title('red channel')
xlabel('Frame t [ms]');
ylabel('Depth (image line time)');
zlabel('Mean intensity at various depths I [%]');

subplot(1,2,2)
surf(X,Y,lineFramesHblu,'lineStyle','none');
colorbar;
title('blue channel')
xlabel('Frame t [ms]');
ylabel('Depth (image line time)');
zlabel('Mean intensity at various depths I [%]');
sgtitle('Comparison of channel mean intensities for all frames over depth');