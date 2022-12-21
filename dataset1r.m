% AMI Project 2022 - Rebecca
% Dataset 1 - H-scan
%--------------------------------------------------------------------------
% load dataset 1 only image contraction at full activation vs full rest
% should be conducted on one frame with only the H-scan for both full
% activation and at rest.

% See separate script dataset1r.m
% --------
% AT REST
% --------
load('C:\Users\Rebecca Viklund\Desktop\AMI project\AMI\Dataset_1\CG_rest_2_omg2.mat');
frame=floor(30/2); % Chose a frame, for dataset 1 one is sufficient for the report.
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

% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

%% Plot one spectrum of one imageline for comparison with the GH spectra
line=floor(size(RF,2)/2); % Image line = 64 out or 128

[pxxlo,flo] = pwelch(GHlo, hamming(512));
f_VECTlo = linspace(0,Fs/2,length(pxxlo));
p_NORMlo = 0.5*pxxlo./max(pxxlo);

[pxxhi,fhi] = pwelch(GHhi,hamming(512));
f_VECThi = linspace(0,Fs/2,length(pxxhi));
p_NORMhi = 0.5*pxxhi./max(pxxhi);

imLineRF = double(squeeze(RF(:,line))); 
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
ylim([	min([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25...
		max([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25]);

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

% rb colorcoding
BmodesrgbHr = zeros(shape(1),shape(2),3);
Bmodesrgb = zeros(shape(1),shape(2),3);
Bmodesrgblo = zeros(shape(1),shape(2),3);
Bmodesrgbhi = zeros(shape(1),shape(2),3);

BmodesrgbHr(:,:,1) = (Bmodeshi(:,:)-min(Bmodeshi(:,:),[],'all')) / (max(Bmodeshi(:,:),[],'all')-min(Bmodeshi(:,:),[],'all')); %Bmodeshi(:,:)./mean2(Bmodeshi(:,:));
BmodesrgbHr(:,:,3) = (Bmodeslo(:,:)-min(Bmodeslo(:,:),[],'all')) / (max(Bmodeslo(:,:),[],'all')-min(Bmodeslo(:,:),[],'all'));
Bmodesrgb(:,:,2) = (Bmodes(:,:)-min(Bmodes(:,:),[],'all')) / (max(Bmodes(:,:),[],'all')-min(Bmodes(:,:),[],'all'));

Bmodesrgbhi(:,:,1) = BmodesrgbHr(:,:,1);
Bmodesrgblo(:,:,3) = BmodesrgbHr(:,:,3);

BmodesrgbHr(:,:,1) = medfilt2(BmodesrgbHr(:,:,1));
BmodesrgbHr(:,:,3) = medfilt2(BmodesrgbHr(:,:,3));
Bmodesrgb = medfilt2(Bmodesrgb(:,:,2));
Bmodesrgbhi = medfilt2(Bmodesrgbhi(:,:,1));
Bmodesrgblo = medfilt2(Bmodesrgblo(:,:,3));

%% analyze intensity in a line; over signal time (depth)
nodepths=size(BmodesrgbHr,1);
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
plot(depths,lineHred,'r');
plot(depths,lineHblu,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channel intensities for frame ',num2str(frame),' over line ',num2str(line)]);
legend({'Red channel','Blue channel'});
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'Red channel','Blue channel','Red trend','Blue trend'});

%% analyze intensity in one frame; over signal time (depth)
nolines=size(BmodesrgbHr,2);
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

%% plot H-scan filtering results
draw_pic(Bmodes(:,:), Bmodeslo, Bmodeshi, [] , 0.05,4);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbHr,0.05,5);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes,BmodesrgbHr,0.05,6);

%% --------
% AT WORK
% --------
load('C:\Users\Rebecca Viklund\Desktop\AMI project\AMI\Dataset_1\CG_contraction_1_omg2.mat');
frame=floor(30/2); % Chose a frame, for dataset 1 one is sufficient for the report.
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

% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

%% Plot one spectrum of one imageline for comparison with the GH spectra
line=floor(size(RF,2)/2); % Image line = 64 out or 128

[pxxlo,flo] = pwelch(GHlo, hamming(512));
f_VECTlo = linspace(0,Fs/2,length(pxxlo));
p_NORMlo = 0.5*pxxlo./max(pxxlo);

[pxxhi,fhi] = pwelch(GHhi,hamming(512));
f_VECThi = linspace(0,Fs/2,length(pxxhi));
p_NORMhi = 0.5*pxxhi./max(pxxhi);

imLineRF = double(squeeze(RF(:,line))); 
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
ylim([	min([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25...
		max([p_NORMlo p_NORMhi p_NORM],[],'all')*1.25]);

%% ---------------------------
% 2D H-scan conv and rgb encoding
%-----------------------------
nolines=size(RF,2);
% convolution
for k=1:nolines
	RFlo(:,k)=conv(RF(:,k),GHlo,'same')./sqrt(Elo);
	RFhi(:,k)=conv(RF(:,k),GHhi,'same')./sqrt(Ehi);
end

% envelope detection
Bmodeslo=sqrt(abs(hilbert(RFlo)));
Bmodeshi=sqrt(abs(hilbert(RFhi)));

% rb colorcoding
BmodesrgbH = zeros(shape(1),shape(2),3);
Bmodesrgb = zeros(shape(1),shape(2),3);
Bmodesrgblo = zeros(shape(1),shape(2),3);
Bmodesrgbhi = zeros(shape(1),shape(2),3);

BmodesrgbH(:,:,1) = (Bmodeshi(:,:)-min(Bmodeshi(:,:),[],'all')) / (max(Bmodeshi(:,:),[],'all')-min(Bmodeshi(:,:),[],'all')); %Bmodeshi(:,:)./mean2(Bmodeshi(:,:));
BmodesrgbH(:,:,3) = (Bmodeslo(:,:)-min(Bmodeslo(:,:),[],'all')) / (max(Bmodeslo(:,:),[],'all')-min(Bmodeslo(:,:),[],'all'));
Bmodesrgb(:,:,2) = (Bmodes(:,:)-min(Bmodes(:,:),[],'all')) / (max(Bmodes(:,:),[],'all')-min(Bmodes(:,:),[],'all'));

Bmodesrgbhi(:,:,1,:) = BmodesrgbH(:,:,1,:);
Bmodesrgblo(:,:,3,:) = BmodesrgbH(:,:,3,:);

% BmodesrgbH(:,:,1) = medfilt2(BmodesrgbH(:,:,1));
% BmodesrgbH(:,:,3) = medfilt2(BmodesrgbH(:,:,3));
% Bmodesrgb = medfilt2(Bmodesrgb(:,:,2));
% Bmodesrgbhi = medfilt2(Bmodesrgbhi(:,:,1));
% Bmodesrgblo = medfilt2(Bmodesrgblo(:,:,3));

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
lines=1:nolines;

for k=1:nodepths
    for j=1:nolines
        Tempr(j)=Bmodesrgbhi(k,j,1);
        Tempb(j)=Bmodesrgblo(k,j,3);
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

%% plot H-scan filtering results
draw_pic(Bmodes(:,:), Bmodeslo, Bmodeshi, [] , 0.05,4);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbH,0.05,5);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes,BmodesrgbH,0.05,6);

%% ---------------
% 1c Comparing muscle at rest and at work H-scan images
%-----------------
figure
subplot(1,2,1)
histogram(BmodesrgbHr,20)
title('At rest')
ylabel('Culmutative sum of occurence');
xlabel('Normalized intensities I [1]')
subplot(1,2,2)
histogram(BmodesrgbH,20)
title('Contraction')
ylabel('Culmutative sum of occurence');
xlabel('Normalized intensities I [1]')
sgtitle('H-scan intensity content')

%%
draw_pic2(BmodesrgbHr,BmodesrgbH,0.05,11);
