% H-scan
% TODO: 
% * Check H2 and H8 polynoms if correct
% * Selection of b1 and b2:
%   Other parameters that should be part of the GH? tau, a and b? Check the
%   papers that described the technique
%   The filters seem to be very broad compared to figures in the papers, 
%   can they be made more narrow?
% * Is filtering needed for the final results? Now rather spike:y. Maybe
%   medfilt is needed or conv2?


load('U:\Till AMI\CG_rest_2_omg2.mat');

frame = 10; % V�lj en frame. F�r dataset 1 beh�ver ni bara en frame. F�r de andra dataset:en kan ni k�ra Hscan p� varje frame, dvs bara loop:a �ver olika frames.
RF = double(RF(1:2177,:,frame)); % H�rdkodat 1:2177 som motsvarar 4cm m�tdjup i detta data.

% TGC - Detta beh�vs bara p� Dataset 1. TGC = "time gain compensation"
% Den f�rst�rker signalen proportionellt mot m�tdjupet pga ljudenergi
% f�rsvinner med �kat m�tdjup...
TGC_VECT = linspace(1,10,size(RF,1))';
TGC_MAT = repmat(TGC_VECT,[1 size(RF,2)]);

RF = RF.*TGC_MAT;

% H�r k�rs en "beamforming". Detta beh�vs ocks� bara p� dataset 1.
% Beh�vs pga data samlas in med "plane-wave imaging" och beamforming skapar
% syntetiskt fokusering.

BF = beamformDAQ(RF,128);
BF = BF - repmat(mean(BF,1),[size(BF,1) 1]);

% Use the beamformed version of RF for all computations below
RF = BF;


Fs = 35; %MHz Samplingsfrekvens av ultraljudet (inte "frame rate"/bilder per sekund)
% Bildens dimensioner motsvarar 
% 40 cm x 40 cm (djup x bredd) f�rdelat p� 2177 x 128 pixlar i RF bilden

%%

% ********************************************************
% Definiera GH2 och GH8
% ********************************************************
% Utmaning kvar = normalisera med energierna ("E") f�r varje GH 
% och optimera b1 och b2

T_duration = 10; % microseconds (the time intervall for the GH pulses)
t = linspace(-T_duration,T_duration,2*T_duration*Fs);
b1 = 0.07;
H2 = (4.*(t./b1).^2 - 2);
GH2 = exp(-(t./(b1)).^2).*H2;
GH2 = GH2./sum(GH2(:));

% Jag testade detta som en variant till GH2 (som f�r mer lokaliserad
% frekvensband). Dvs: Att anv�nda GH2 = GH8 men med olika b-skalning s� att
% det blir en h�g resp l�gbandpass.
b2 = 0.12;
H8 = (256.*(t./b2).^8 - 3585.*(t./b2).^6 + 13440.*(t./b2).^4 - 13440.*(t./b2).^2 + 1680);
GH8 = exp(-(t./(b2)).^2).*H8;
GH8 = GH8./sum(GH8(:));
GH2 = GH8;

[pxx2,f2] = pwelch(GH2, hamming(512));
f_VECT2 = linspace(0,Fs/2,length(pxx2));

b2 = 0.07;
H8 = (256.*(t./b2).^8 - 3585.*(t./b2).^6 + 13440.*(t./b2).^4 - 13440.*(t./b2).^2 + 1680);
GH8 = exp(-(t./(b2)).^2).*H8;
GH8 = GH8./sum(GH8(:));
[pxx8,f1] = pwelch(GH8,hamming(512));
f_VECT8 = linspace(0,Fs/2,length(pxx8));


% *******************************************
% Plot examples
% *******************************************
figure(2);
subplot(2,1,1); plot(t, GH2, 'r'); ylabel('GH_2(t)');
subplot(2,1,2); plot(t, GH8, 'b'); ylabel('GH_8(t)'); xlabel('Time, us');

figure(1); clf; hold on;
plot(f_VECT2, 0.5*pxx2./max(pxx2), 'r', 'LineWidth', 2);
plot(f_VECT8, 0.5*pxx8./max(pxx8), 'b', 'LineWidth', 2);
xlabel('Frequency, MHz');
ylabel('Amplitude, n.u.');

% Plot one spectrum of one imageline for comparison with the GH spectra
imLineRF = double(squeeze(RF(:,64))); % Image line = 64 out or 128
[pxx,f] = pwelch(imLineRF);
f_VECT = linspace(0,Fs/2,length(f));

figure(1);
plot(f_VECT, sqrt(pxx)./max(sqrt(pxx)), 'k', 'LineWidth', 2);
legend({'GH_2','GH_8','RF imageline'});

%%

% H�r gissar jag bara, ni f�r kolla hur det ska ber�knas

% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(imLineRF)));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
H2 = conv(imLineRF, GH2, 'same'); % Obs: 'same' argument is important
% envelope conversion (I was guessing that this is needed, please double-check)
H2 = sqrt(abs(hilbert(H2)));

H8 = conv(imLineRF, GH8, 'same');
H8 = sqrt(abs(hilbert(H8)));

Red_ImLine = H2./H8; % St�mmer det? Jag tror dessa varianter anv�ndes i n�gon av artiklarna
Blue_ImLine = H8./H2; % St�mmer det?

figure(11); clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
legend({'Red channel','Blue channel'});
