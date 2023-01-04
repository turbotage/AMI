% AMI Project 2022 - Rebecca
% Dataset 1 - H-scan
%--------------------------------------------------------------------------
% load dataset 1 only image contraction at full activation vs full rest
% should be conducted on one frame with only the H-scan for both full
% activation and at rest.

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
TGC_V = linspace(1,10,size(RF,1))';
TGC_M = repmat(TGC_V,[1 size(RF,2)]);

RF = RF.*TGC_M;

% Here a "beamforming" is conducted. This is only needed for dataset 1.
% This is due to the collection method "plane-wave imaging" used which
% requires synthetic focusing which is provided by beamforming techniques.

BF = beamformDAQ(RF,128);
BF = BF - repmat(mean(BF,1),[size(BF,1) 1]);

% Use the beamformed version of RF for all computations below
RF = BF;

RF = single(RF(:,:));

%% signal profile
%  analyzes mean intensity off all frames; over signal time (depth)

for k=1:size(RF,1)
    for j=1:size(RF,2)
        Temp(j)=RF(k,j);
    end
    line(k)=mean(Temp);
end

line(isnan(line))=0;

figure
plot(line)
title('Dataset 1 CG rest mean signal profile')
xlabel('Depth (image line time)');
ylabel('Mean intensity at various depths I [%]');

%%
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

% GH high pass
b2 = 0.125;
ordhi = 32;
Hhi = hermiteH(ordhi, t./b2); % order 32
GHhi = exp(-(t./(b2)).^2).*Hhi;

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

% computing analytic energies for conv
Elo=prod(1:2:(2*ordlo-1))*sqrt(pi/2);
Ehi=prod(1:2:(2*ordhi-1))*sqrt(pi/2);

%% Plot one spectrum of one imageline for comparison with the GH spectra
line=floor(size(RF,2)/2); % Image line = 64 out or 128

[pxxlo,flo] = pwelch(GHlo, hamming(512));
f_VECTlo = linspace(0,Fs/2,length(flo));
p_NORMlo = 0.5*pxxlo./max(pxxlo);

[pxxhi,fhi] = pwelch(GHhi,hamming(512));
f_VECThi = linspace(0,Fs/2,length(fhi));
p_NORMhi = 0.5*pxxhi./max(pxxhi);

imLineRF = double(squeeze(RF(:,line))); 
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));
p_NORM = sqrt(pxx)./max(sqrt(pxx));

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

BmodesrgbHr(:,:,1) = myrgbencoder(Bmodeslo);
BmodesrgbHr(:,:,3) = myrgbencoder(Bmodeshi);
Bmodesrgb(:,:,2) = myrgbencoder(Bmodes);

Bmodesrgblo(:,:,1) = BmodesrgbHr(:,:,1);
Bmodesrgbhi(:,:,3) = BmodesrgbHr(:,:,3);


BmodesrgbHr(:,:,1) = medfilt2(BmodesrgbHr(:,:,1),[3,22]);
BmodesrgbHr(:,:,3) = medfilt2(BmodesrgbHr(:,:,3),[3,22]);
Bmodesrgb(:,:,2) = medfilt2(Bmodesrgb(:,:,2),[3,22]);

Bmodesrgblo(:,:,1) = BmodesrgbHr(:,:,1);
Bmodesrgbhi(:,:,3) = BmodesrgbHr(:,:,3);

%% analyze intensity in a line; over signal time (depth)
nodepths=size(BmodesrgbHr,1);
depths=1:nodepths;

for j=1:nodepths
    lineHlo(j)=Bmodesrgblo(j,line,1);
    lineHhi(j)=Bmodesrgbhi(j,line,3);
end

lineHlo(isnan(lineHlo))=0;
lineHhi(isnan(lineHhi))=0;

lineHlop=lineHlo./lineHhi;
lineHhip=lineHhi./lineHlo;

figure(3); clf; hold on;
plot(depths,lineHlop,'r');
plot(depths,lineHhip,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of filter channel intensities for frame\newline',num2str(frame),' over line ',num2str(line)]);
legend({'Red channel','Blue channel'});
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'High pass','Low pass'});

%% analyze intensity in one frame; over signal time (depth)
nolines=size(BmodesrgbHr,2);
lines=1:nolines;

for k=1:nodepths
    for j=1:nolines
        Templo(j)=Bmodesrgblo(k,j,1);
        Temphi(j)=Bmodesrgbhi(k,j,3);
    end
    lineFrameHlo(k)=mean(Templo);
    lineFrameHhi(k)=mean(Temphi);
end

lineFrameHlo(isnan(lineFrameHlo))=0;
lineFrameHhi(isnan(lineFrameHhi))=0;

lineFrameHlop=lineFrameHlo./lineFrameHhi;
lineFrameHhip=lineFrameHhi./lineFrameHlo;

figure(3); clf; hold on;
plot(depths,lineFrameHhip,'r');
plot(depths,lineFrameHlop,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channel mean intensities for\newlineframe ',num2str(frame),' over depth']);
legend({'Red channel','Blue channel'});
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'Red channel','Blue channel'});

%% manual climc selection
figure
histogram(BmodesrgbHr,20)
title('At rest intensity distribution')
ylabel('Culmutative sum of occurence');
xlabel('Normalized intensities I [1]')

climcr=[];

%% plot H-scan filtering results
draw_pic(Bmodes,Bmodeslo,Bmodeshi, [], 0.05, []);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbHr, 0.05, climcr);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes,BmodesrgbHr, 0.05, climcr);

%% --------
% AT WORK
% --------
load('C:\Users\Rebecca Viklund\Desktop\AMI project\AMI\Dataset_1\CG_contraction_1_omg2.mat');
frame=floor(30/2); % Chose a frame, for dataset 1 one is sufficient for the report.
% ----------------------
% Pre-processing of data
%-----------------------
RF = double(RF(1:2177,:,frame)); %1:2177 corresp to 4cm depth in this dataset..

% TGC - "time gain compensation"
% It amplifies the signal proportional to the US depth since it weakens as
% it interacts with tissue.
TGC_V = linspace(1,10,size(RF,1))';
TGC_M = repmat(TGC_V,[1 size(RF,2)]);

RF = RF.*TGC_M;

% Here a "beamforming" is conducted. This is only needed for dataset 1.
% This is due to the collection method "plane-wave imaging" used which
% requires synthetic focusing which is provided by beamforming techniques.

BF = beamformDAQ(RF,128);
BF = BF - repmat(mean(BF,1),[size(BF,1) 1]);

% Use the beamformed version of RF for all computations below
RF = BF;

RF= single(squeeze(RF(:,:,:)));

%% signal profile
%  analyzes mean intensity off all frames; over signal time (depth)

for k=1:size(RF,1)
    for j=1:size(RF,2)
        Temp(j)=RF(k,j);
    end
    line(k)=mean(Temp);
end

line(isnan(line))=0;

figure
plot(line,'-r')
title('Dataset 1 CG contraction mean signal profile')
xlabel('Depth (image line time)');
ylabel('Mean intensity at various depths I [%]');

%%
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

% GH high pass
b2 = 0.135;
ordhi = 32;
Hhi = hermiteH(ordhi, t./b2); % order 32
GHhi = exp(-(t./(b2)).^2).*Hhi;

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
plot(f_VECTlo, p_NORMlo, 'r', 'LineWidth', 1.5);
plot(f_VECThi, p_NORMhi, 'b', 'LineWidth', 1.5);
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
BmodesrgbH = zeros(shape(1),shape(2),3);
Bmodesrgb = zeros(shape(1),shape(2),3);
Bmodesrgblo = zeros(shape(1),shape(2),3);
Bmodesrgbhi = zeros(shape(1),shape(2),3);


BmodesrgbH(:,:,1) = myrgbencoder(Bmodeslo);
BmodesrgbH(:,:,3) = myrgbencoder(Bmodeshi);
Bmodesrgb(:,:,2) = myrgbencoder(Bmodes);

Bmodesrgblo(:,:,1) = BmodesrgbH(:,:,1);
Bmodesrgbhi(:,:,3) = BmodesrgbH(:,:,3);


BmodesrgbH(:,:,1) = medfilt2(BmodesrgbH(:,:,1),[3,22]);
BmodesrgbH(:,:,3) = medfilt2(BmodesrgbH(:,:,3),[3,22]);
Bmodesrgb(:,:,2) = medfilt2(Bmodesrgb(:,:,2),[3,22]);

Bmodesrgblo(:,:,1) = BmodesrgbH(:,:,1);
Bmodesrgbhi(:,:,3) = BmodesrgbH(:,:,3);
%% analyze intensity in a line; over signal time (depth)
nodepths=size(BmodesrgbH,1);
depths=1:nodepths;

for j=1:nodepths
    lineHlo(j)=Bmodesrgblo(j,line,1);
    lineHhi(j)=Bmodesrgbhi(j,line,3);
end

lineHlo(isnan(lineHlo))=0;
lineHhi(isnan(lineHhi))=0;

lineHlop=lineHlo./lineHhi;
lineHhip=lineHhi./lineHlo;

figure(3); clf; hold on;
plot(depths,lineHlop,'r');
plot(depths,lineHhip,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of filter channel intensities for frame\newline',num2str(frame),' over line ',num2str(line)]);
legend({'Red channel','Blue channel'});
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'High pass','Low pass'});

%% analyze intensity in one frame; over signal time (depth)
nolines=size(BmodesrgbH,2);
lines=1:nolines;

for k=1:nodepths
    for j=1:nolines
        Templo(j)=Bmodesrgblo(k,j,1);
        Temphi(j)=Bmodesrgbhi(k,j,3);
    end
    lineFrameHlo(k)=mean(Templo);
    lineFrameHhi(k)=mean(Temphi);
end

lineFrameHlo(isnan(lineFrameHlo))=0;
lineFrameHhi(isnan(lineFrameHhi))=0;

lineFrameHlop=lineFrameHlo./lineFrameHhi;
lineFrameHhip=lineFrameHhi./lineFrameHlo;

figure(3); clf; hold on;
plot(depths,lineFrameHlop,'r');
plot(depths,lineFrameHhip,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of filtered channel mean intensities for\newlineframe ',num2str(frame),' over depth']);
xlabel('Depth (image line time)'); ylabel('Channel intensity I [%]');
legend({'Low pass','High pass'});

%% manual climc selection
figure
histogram(BmodesrgbH,20)
title('Contracted intensity distribution')
ylabel('Culmutative sum of occurence');
xlabel('Normalized intensities I [1]')

climc=[];

%% plot H-scan filtering results
draw_pic(Bmodes,Bmodeslo, Bmodeshi, [], 0.05, []);

%% plot 2D colorcoded H-scan filtering results
draw_pic(Bmodesrgb,Bmodesrgblo,Bmodesrgbhi,BmodesrgbH, 0.05, climc);

%% plot comparison between B-mode and H-scan
draw_pic2(Bmodes,BmodesrgbH, 0.05, climc);

%% ---------------
% 1c Comparing muscle at rest and at work H-scan images
%-----------------
climtog=[min([climr clim],[],'all'),...
		max([climr clim],[],'all')];
draw_pic2(BmodesrgbHr,BmodesrgbH, 0.05, climctog);
