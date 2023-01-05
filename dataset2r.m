% AMI Project 2022 - Rebecca
% Dataset 2 - H-scan
%--------------------------------------------------------------------------
% load dataset 2 stimulated contraction w ~1.6Hz i.e. 2-3 Hz on 1 muscle complex
% should be conducted on several frames and is to be filtered using
% svd and H-scan.

% ----------------------
% Pre-processing data
%-----------------------
% Specify dataset
chosen_dataset='0913';
if strcmp(chosen_dataset,'1555')
    stimFreq = '~8Hz';

elseif strcmp(chosen_dataset,'1311')
    stimFreq= '~7.5Hz';


elseif strcmp(chosen_dataset,'0913')
    stimFreq= '2-3Hz';

else
    disp('That dataset does not exist')
end

%% read data
dspath=strcat('C:\Users\Rebecca Viklund\Desktop\AMI project\AMI\Dataset_2\181023_',chosen_dataset);

if strcmp(chosen_dataset,'1311')
    dspath1=strcat(dspath,'_rs.mat');
    rfmat = single(load(dspath1).rfmat_downsampled);
	b1=.16;
	b2=.135;

else
    dspath1=strcat(dspath,'_rsf.mat');
    rfmat = single(load(dspath1).rfdat);
	b1=.145;
	b2=.125;
end

rfmat = rfmat(1:1000,:,:);
tgc_vect = linspace(1,5,size(rfmat,1))';
tgc_mat = repmat(tgc_vect,[1 size(rfmat,2)]);
rfmat=rfmat.*tgc_mat;

dspath2=strcat(dspath,'_ROI.mat');
roimat = load(dspath2).ROI;

Bmodes = sqrt(abs(hilbert(squeeze(rfmat(:,:,:)))));
shape = size(Bmodes);

%% signal profile
%  analyzes mean intensity off all frames; over signal time (depth)
clear X Y lineFrames
for l=1:size(rfmat,3)
    for k=1:size(rfmat,1)
        for j=1:size(rfmat,2)
            Temp(j)=rfmat(k,j,l);
        end
    lineFrames(k,l)=mean(Temp);
    end
end

lineFrames(isnan(lineFrames))=0;

[X,Y] = meshgrid(1:size(rfmat,3),1:size(rfmat,1));
figure
surf(X,Y,lineFrames,'lineStyle','none');
title(['Dataset 2 ', chosen_dataset, ' ', stimFreq, '\newlinestimulation mean signal profile'])
xlabel('Frame t [0.5ms]');
ylabel('Depth (image line time)');
zlabel('Mean intensity at various depths I [%]');
colorbar;

%% ------------------------
% SVD vs TVI - Viktors part
% -------------------------
% -----------------------
% Pre-processing TVI data 
%------------------------
% fn_STR = '181023_1311tvi_rs.mat';
% h = matfile(fn_STR);
% tvif_d = single(h.tvi_downsampled(:,:,1:250));
% tvif_d = sqrt(abs(hilbert(squeeze(tvif_d(:,:,:)))));
%------------------------
% Filtering and comparing
%------------------------
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
%----------------------
% checking filter settings in 1D by taking a line/col of the img for
% freq analysis
Fs = 40; % MHz Samplingfreq of the US not the signal
T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);

% GH low pass
ordlo = 2;
Hlo = hermiteH(ordlo, t./b1);
GHlo = exp(-(t./(b1)).^2).*Hlo;

% GH high pass
ordhi = 16;
Hhi = hermiteH(ordhi, t./b2); % order 32
GHhi = exp(-(t./(b2)).^2).*Hhi;

%%
figure(1);
subplot(2,1,1)
plot(t,GHlo,'r', 'LineWidth', 1.5);
title(['GH_{',num2str(ordlo),'}'])
ylim([min(GHlo,[],'all')*1.25 max(GHlo,[],'all')*1.25]);
subplot(2,1,2)
plot(t,GHhi,'b', 'LineWidth', 1.5);
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

%%
figure(2); clf; hold on;
plot(f_VECT, p_NORM,'-','color',[0 .5 0],'LineWidth', 1.5);
plot(f_VECTlo, p_NORMlo, 'r', 'LineWidth', 1.5);
plot(f_VECThi, p_NORMhi, 'b', 'LineWidth', 1.5);
sgtitle(['PSD for frame ',num2str(frame),' and line ',num2str(line)]);
xlabel('Frequency f [MHz]'); ylabel('Amplitude A [1]');
legend({'RF imageline','Low pass','High pass'});
ylim([	min([p_NORMlo p_NORMhi [p_NORM;ones(length(p_NORMlo)-length(p_NORM),1)]],[],'all')*1.15...
		max([p_NORMlo p_NORMhi [p_NORM;zeros(length(p_NORMlo)-length(p_NORM),1)]],[],'all')*1.15]);

%% ---------------------------
% 2D H-scan conv and rgb encoding
%-----------------------------
clear RF_MATlo RF_MAThi
clear Bmodeslo Bmodeshi
% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

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

%% rb colorcoding and filtering
BmodesrgbH = zeros(shape(1),shape(2),3,noframes);
Bmodesrgb = zeros(shape(1),shape(2),3,noframes);
Bmodesrgblo = zeros(shape(1),shape(2),3,noframes);
Bmodesrgbhi = zeros(shape(1),shape(2),3,noframes);

BmodesrgbH(:,:,1,:) = myrgbencoder(Bmodeslo);
BmodesrgbH(:,:,3,:) = myrgbencoder(Bmodeshi);
Bmodesrgb(:,:,2,:) = myrgbencoder(Bmodes);

BmodesrgbH(:,:,1,:) = mymedfilt(BmodesrgbH(:,:,1,:), [10,3]);
disp('1')
BmodesrgbH(:,:,3,:) = mymedfilt(BmodesrgbH(:,:,3,:), [10,3]);
disp('2')
Bmodesrgb(:,:,2,:) = mymedfilt(Bmodesrgb(:,:,2,:), [10,3]);
disp('3')

Bmodesrgblo(:,:,1,:)=BmodesrgbH(:,:,1,:);
Bmodesrgbhi(:,:,3,:)=BmodesrgbH(:,:,3,:);

%% manual climc selection
for j=1:noframes
    rgblomean(j)=mean2(Bmodesrgblo(:,:,:,j));
    rgbhimean(j)=mean2(Bmodesrgbhi(:,:,:,j));
end

t = tiledlayout(1,2,'TileSpacing','Compact');

ax1=nexttile;
histogram(rgblomean,20,'FaceColor','red');
%ax1.XTick=[.0382 .03934 .03984 .040496 .041316];
title(ax1,'Low pass');

ax2=nexttile;
histogram(rgbhimean,20,'FaceColor','blue');
%ax2.XTick=[.0453 .045965 .046497 .047162 0.047827];
title(ax2,'High pass');

ylabel(t,'Culmutative sum of occurence');
xlabel(t,'Normalized intensities I [1]');
title(t,'H-scan mean channel intensity content for dataset 3')

climc=[0.03,0.04];

%% plot H-scan filtering results
draw_pic(Bmodes, Bmodeslo, Bmodeshi, [], 0.05, []);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbH, 0.05, climc);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes,BmodesrgbH, 0.1, climc);

%% analyze intensity in a pixel; over measurement time (frames)
depth=floor(size(BmodesrgbH,1)/2);

for j=1:noframes
    pixelHlo(j)=Bmodesrgblo(depth,line,1,j);
    pixelHhi(j)=Bmodesrgbhi(depth,line,3,j);
end

pixelHhi(isnan(pixelHhi))=0;
pixelHlo(isnan(pixelHlo))=0;

pixelHlo=pixelHlo./pixelHhi;
pixelHhi=pixelHhi./pixelHlo;

[pks,locs] = findpeaks(pixelHlo,frames);
[pks1,locs1] = findpeaks(pixelHhi,frames);

pixelHloMax=makima(locs,pks,frames);
pixelHhiMax=makima(locs1,pks1,frames);

figure; clf; hold on;
plot(frames,pixelHlo','-r');
plot(frames,pixelHhi,'-b');
plot(frames,pixelHloMax,'-r','LineWidth',1.5);
plot(frames,pixelHhiMax,'-b','LineWidth',1.5);
%ylim([0 1]);
xlabel('Frames t [ms]'); ylabel('Channel intensity I [%]');
sgtitle(['Comparison of filter channel intensities in pixel (',num2str(line),'x',num2str(depth),')']);
legend({'Low pass','High pass','Red trend','Blue trend'});

%% analyze std of intensity in each frame; over measurement time (frames)

for j=1:noframes
    frameHlo(:,j) = std2(Bmodesrgblo(:,:,1,j));
    frameHhi(:,j) = std2(Bmodesrgbhi(:,:,3,j));
end

frameHhi(isnan(frameHhi))=0;
frameHlo(isnan(frameHlo))=0;

frameHlo=frameHlo./frameHhi;
frameHhi=frameHhi./frameHlo;

figure; clf; hold on;
plot(frames,frameHlo','-r');
plot(frames,frameHhi,'-b');
% ylim([0 1]);
xlabel('Frames t [ms]'); ylabel('Channel intensity std \sigma [1]');
sgtitle(['Comparison of filter channel std of intensities in frames']);
legend({'Low pass','High pass','Red trend','Blue trend'});

%% analyze intensity in a line; over signal time (depth)
nodepths=size(BmodesrgbH,1);
depths=1:nodepths;

for j=1:nodepths
    lineHlo(j)=Bmodesrgblo(j,line,1,frame);
    lineHhi(j)=Bmodesrgbhi(j,line,3,frame);
end

lineHhi(isnan(lineHhi))=0;
lineHlo(isnan(lineHlo))=0;

lineHhip=lineHhi./lineHlo;
lineHlop=lineHlo./lineHhi;

figure(3); clf; hold on;
plot(depths,lineHlop,'r');
plot(depths,lineHhip,'b');
sgtitle(['Comparison of filter channel intensities for frame ',num2str(frame),' over line ',num2str(line)]);
legend({'Red channel','Blue channel'});
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'Low pass','High pass','Red trend','Blue trend'});

%% analyze intensity in one frame; over signal time (depth)
nolines=size(BmodesrgbH,2);
lines=1:nolines;

for k=1:nodepths
    for j=1:nolines
        Temphi(j)=Bmodesrgbhi(k,j,3,frame);
        Templo(j)=Bmodesrgblo(k,j,1,frame);
    end
    lineFrameHhi(k)=mean(Temphi);
    lineFrameHlo(k)=mean(Templo);
end

lineFrameHhi(isnan(lineFrameHhi))=0;
lineFrameHlo(isnan(lineFrameHlo))=0;

lineFrameHhip=lineFrameHhi./lineFrameHlo;
lineFrameHlop=lineFrameHlo./lineFrameHhi;

figure(3); clf; hold on;
plot(depths,lineFrameHlop,'r');
plot(depths,lineFrameHhip,'b');
sgtitle(['Comparison of channel mean intensities for\newlineframe ',num2str(frame),' over depth']);
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'Low pass','High pass','Red trend','Blue trend'});

%% analyze mean intensity off all frames; over signal time (depth)
clear Templo Temphi lineFramesHlo lineFramesHhi
nodepths=size(BmodesrgbH,1);
depths=1:nodepths;
nolines=size(BmodesrgbH,2);
lines=1:nolines;

for m=1:noframes
    for k=1:nodepths
        for j=1:nolines
            Templo(j)=Bmodesrgblo(k,j,1,m);
            Temphi(j)=Bmodesrgbhi(k,j,3,m);
        end
    lineFramesHhi(k,m)=mean(Temphi);
    lineFramesHlo(k,m)=mean(Templo);
    end
end

lineFramesHhi(isnan(lineFramesHhi))=0;
lineFramesHlo(isnan(lineFramesHlo))=0;

lineFramesHlop=lineFramesHlo./lineFramesHhi;
lineFramesHhip=lineFramesHhi./lineFramesHlo;

[X,Y] = meshgrid(frames,depths);
figure
g = gcf;
g.WindowState = 'maximized';

t = tiledlayout(1,2,'TileSpacing','Compact');
ax3=nexttile;
surf(X,Y,lineFramesHlop,'lineStyle','none');
title('Low pass channel')
xlabel('Frame t [ms]');
ylabel('Depth (image line time)');

ax3=nexttile;
surf(X,Y,lineFramesHhip,'lineStyle','none');
title('High pass channel')
xlabel('Frame t [ms]');
ylabel('Depth (image line time)');

% 14x9 cm export setup
ylabel(t,'Mean intensity at various depths I [%]','FontSize',8);

title(t,['H-scan channel mean intensities for dataset 2\newline', chosen_dataset, ' ', stimFreq, ' stimulation signal']);
cb = colorbar;
cb.Layout.Tile = 'east';

%% Mask the image
% for j=1:noframes
%     Hmasked(:,:,:,j) = BmodesrgbH(:,:,:,j).*repmat(roimat,[1,1,3]);
% end
% 
% % plot mask vs H-scan
% opacity=.1;
% figure
% W1=BmodesrgbH(:,:,:,1);
% W2=roimat;
% W3=Hmasked(:,:,:,1);
% 
% subplot(1,2,1)
% img1=image(W1,'CDataMapping','scaled');
% hold on;
% img2=image(W2,'CDataMapping','scaled','AlphaData',opacity*ones(size(roimat)));colormap gray;
% title(['H-scan with ROI overlapped at opacity ',num2str(opacity)])
% subplot(1,2,2)
% img3=image(W3,'CDataMapping','scaled');
% title('masked out image area');
% pause(2)
% for j=1:floor(noframes/4)
%     W1 = BmodesrgbH(:,:,:,j);
%     W1 = imadjust(W1, [0,1],[]);
% 
%     W3 = Hmasked(:,:,:,j);
%     W3 = imadjust(W3, [0,1],[]);
%     
%     set(img1, 'CData', W1);
%     set(img2, 'CData', W2);
%     set(img3, 'CData', W3);
% 
%     %drawnow limitrate;
%     drawnow();
%     clim([0,1])
%     pause(.05);
%     disp(j);
% end