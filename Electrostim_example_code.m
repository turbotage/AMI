% Test E-stim dataset #2


% Luke - bra/dålig ?????
% h = matfile('U:\Till AMI\181023_1526.mat'); 

% CG, finns något mitt i bild, men inte steady state i stimulering
% h = matfile('V:\181023_1556.mat'); 
% CG, finns något mitt i bild, men inte steady state i stimulering


% Tomas går ej se något tydligt
% h = matfile('V:\181023_1013.mat'); 

% Tomas går ej se något tydligt
% h = matfile('V:\181023_1013.mat'); 

% h = matfile('V:\181023_1013.mat'); 

% **************************
% USE the ones below!
% **************************

% Karro #2 - ROI OK
% fn_STR = 'V:\181023_0913.mat';
% h = matfile(fn_STR); 
% stim frekvens 2-3 Hz. ganska stabil

% Madeleine #1 - ROI OK
fn_STR = 'C:\Users\chgr05\OneDrive - Umeå universitet\Applied Medical Imaging 2022-P2\181023_1311.mat';
h = matfile(fn_STR); 
% Stim frekvens = 7.5 Hz stabil
% Går se rätt tydligt uppe i högra hörnet och tydligt återkommande

% CG #1 - ROI OK
% fn_STR = 'V:\181023_1555.mat';
% h = matfile(fn_STR); 
% stim frekvens ~8Hz. Stabil

% Load down-sampled version of data (2 sek but at 500Hz rate, only first
% 1000 samples in y, means 2 cm depth, this is enough, skip the rest).
% Does not seem to present aliasing effect, and should not because the
% signal content is <<100Hz.
RF = h.rf_data_set(1:1000,:,1:4:4000);

% Apply TGC
TGC_VECT = linspace(1,10,size(RF,1))';
TGC_MAT = repmat(TGC_VECT,[1 size(RF,2)]);

%%
% figure(1); clf;

% Reconstruct B-mode signal
B = zeros(size(RF));
for n = 1:size(RF,3),
    n
    RFframe = squeeze(RF(:,:,n)).*TGC_MAT;
    B(:,:,n) = sqrt(abs(hilbert(RFframe)));
    
%     subplot(121);
%     imagesc(sign(RFframe).*sqrt(abs(RFframe)));
%     colormap gray
%     subplot(122);
%     imagesc(sqrt(abs(hilbert(RFframe))));
    
%     pause
end

%% Band-pass filtering option
% Bf = zeros(size(B));
% [b,a] = butter(4,[4 30]./(500/2), 'bandpass');
% for r = 1:size(B,1), 
%     r
%     for c = 1:size(B,2),
%         Bf(r,c,:) = filtfilt(b,a,squeeze(B(r,c,:)));
%     end
% end

%% Visualization as sequence
figure(10);
for n = 1:1:size(B,3), 
    TMP = squeeze(B(:,:,n));
    TMP = conv2(TMP,ones(20,2),'same');
    imagesc(TMP(1:1:end,:)); colormap hot; title(n); pause; 
end


%