%  AMI Project 2022

%% load dataset 1 ... only image contraction w ~ Hz full activation vs full rest

%% load dataset 2 stimulated contraction w ~1.6Hz 1 muscle complex
rfmat_dsf = single(load('181023_1311_rs.mat').rfmat_downsampled);

Bmodes = sqrt(abs(hilbert(squeeze(rfmat_dsf(:,:,:)))));
shape = size(Bmodes);
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
clear Bmodes

% svd
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(25:end, 25:end) = 0;
Snew(1:20,1:20) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

% H-scan

% MY COMPUTER CAN'T HANDLE THIS DATASET!

%% load dataset 3 voluntary contraction w ~5-25Hz 2-4 muscle complexes
load("RF_MAT.mat");
RF_MAT= single(squeeze(RF_MAT(:,:,:)));
Bmodes = single(sqrt(abs(hilbert(squeeze(RF_MAT(:,:,:))))));
shape = size(Bmodes);
RF_MAT_f=reshape(RF_MAT, shape(1)*shape(2), shape(3));

%% H-scan
% generate kernel of hermitian polynomials of ord 2 and 8
% with matrices
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

%% RGB code signal and convolute
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

%% Plot rgb signal for 1 line in 1 frame
% Compute "B-mode" image line = Green channel
Green_ImLine = sqrt(abs(hilbert(RF_MAT(:,line,frame))));

% Compute H-scan
% Convolution. This is the "time-domain filtering" procedure using the GH2 and GH8.
H2 = conv(RF_MAT(:,line,frame), GH2, 'same');

% envelope conversion (I was guessing that this is needed, please double-check)
H2 = sqrt(abs(hilbert(H2)));

H8 = conv(RF_MAT(:,line,frame), GH8, 'same');
H8 = sqrt(abs(hilbert(H8)));

Red_ImLine = H8./H2; % Stämmer det? Jag tror dessa varianter användes i någon av artiklarna
Blue_ImLine = H2./H8; % Stämmer det?

figure; clf; hold on;
plot(Red_ImLine,'r');
plot(Blue_ImLine,'b');
xlabel('Depth (image line time)');
sgtitle(['Comparison of channels for frame \newline',num2str(frame),' and line ',num2str(line)]);
legend({'Red channel','Blue channel'});

%% Combination with colorcoding in freq space to RGB and then combined into 
% a new image RGB colormap R for low freqs, B for high, G for original sig
% envelope.

% RGB channeling intensites G=signal, B=GH8, R=GH2
% mymap = [	1 0 0 ;
%			0 1 0 ;
%			0 0 1];

% Colorcode in 2D
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
fn_STR = '181023_1311tvi_rs.mat';
h = matfile(fn_STR);
tvif_d = single(h.tvi_downsampled(:,:,1:250));
tvif_d = sqrt(abs(hilbert(squeeze(tvif_d(:,:,:)))));

%% SVD

Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(25:end, 25:end) = 0;
Snew(1:20,1:20) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%% SVD vs TVI 
% Prefiltering
% median3 and butter order 4 on both TVI and SVD, keeping the freqs in
% range 5-25 Hz
% suggestions
% tripple median filtering
%Bmodes_fin=medfilt3(Bmodes_new);
%tvi_fin=medfilt3(tvif);
% spatial filtering
%C = conv2(A, B8);
% freq filtering
% idea: bandpass on every pixelintensity over activationtime 5-15(30) Hz
% for dataset 2 but around 1.5 Hz for dataset 3
%[b,a] = butter(4,[5 25]./(Fsamp/2),'bandpass');
%y_filt = filtfilt(b,a,y);

[b1,b2] = butter(1, 200/250,'low');
tvif = filter(b1,b2,tvif_d,[],3);
Bmodes_new = filter(b1,b2,Bmodes_new,[],3);

for j=1:noframes
	tviff(:,:,j) = medfilt2(tvi_line_filtered(:,:,j));
	Bmodes_newff(:,:,j) = medfilt2(Bmodes_new(:,:,j));
end

%%
tvi_line_filtered = full_line_filter(tvif);
draw_pic2(tvi_line_filtered, tviff);
%%

after_lf_masks = logical(zeros(size(tvi_line_filtered)));
for i=1:size(tvi_line_filtered,3)
    after_lf_masks(:,:,i) = imbinarize(tvi_line_filtered(:,:,i),"adaptive","ForegroundPolarity","dark");
end

draw_pic2(tvi_line_filtered, after_lf_masks);

%%
tviffconv = zeros(size(tviff));
for i=1:size(tviff,3)
    tviffconv(:,:,i) = conv2(tviff(:,:,i),ones(4,3)/3,'same');
end

draw_pic2(tviff, tviffconv);

%%
stdpic = std(tviffconv,0,3);

zmax = size(tviffconv,3);
steplen = round(zmax/5);
under_step = 1;
sumtvi = zeros(size(tviffconv,1),size(tviffconv,2));
for i=steplen:steplen:zmax
    sumtvi = sumtvi + std(tviffconv(:,:,under_step:i), 0, 3);
    under_step = under_step + steplen;
end

%stdpic = conv2(stdpic,ones(5,5)/25, 'same');

%%
figure;
imagesc(sumtvi);
colormap gray;
axis('square');

figure;
tvisum_mask = imbinarize(sumtvi, "adaptive", "Sensitivity", 0.001);
imagesc(tvisum_mask);
colormap gray;
axis('square');

%%
% Mean intensity of line filtered
tlen = size(tvi_line_filtered, 3);
zmax = size(tvi_line_filtered, 3);
xstart = 1200;
xstop = 1300;
ystart = 80;
ystop = 128;

avg_intense = reshape(mean(mean( ...
    tvi_line_filtered(xstart:xstop,ystart:ystop,:),1),2),[tlen,1]);
figure;
subplot(1,2,1);
plot(avg_intense);

subplot(1,2,2);
fD = fft(avg_intense);
plot(abs(fD));

[safD,idx] = sort(abs(fD),'descend');
%fD(idx(2:15)) = 0;
rmid = 33;
fD((rmid-1):(rmid+1)) = 0; fD((zmax-rmid+1):(zmax-rmid+3)) = 0;
rmid = 17;
fD((rmid-1):(rmid+1)) = 0; fD((zmax-rmid+1):(zmax-rmid+3)) = 0;
rmid = 49;
fD((rmid-1):(rmid+1)) = 0; fD((zmax-rmid+1):(zmax-rmid+3)) = 0;
figure;
subplot(1,2,1);
plot(abs(fD));
ifD = ifft(fD);
subplot(1,2,2);
plot(ifD);

%%
freq_tlf = fft(tviffconv,[],3);
rmid = 33;
freq_tlf(:,:,(rmid-1):(rmid+1)) = 0; freq_tlf(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
rmid = 17;
freq_tlf(:,:,(rmid-1):(rmid+1)) = 0; freq_tlf(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
rmid = 49;
freq_tlf(:,:,(rmid-1):(rmid+1)) = 0; freq_tlf(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;

tlf_pf = ifft(freq_tlf,[],3);

%%
draw_pic2(tvif_d, tlf_pf);

%%
tlf_pf_mean = reshape(mean(mean( ...
    tlf_pf(xstart:xstop,ystart:ystop,:),1),2),[tlen,1]);
figure;
plot(tlf_pf_mean);

%%
zmax = size(tlf_pf,3);
steplen = round(zmax/5);
under_step = 1;
sumtvi = zeros(size(tlf_pf,1),size(tlf_pf,2));
for i=steplen:steplen:zmax
    sumtvi = sumtvi + std(tlf_pf(:,:,under_step:i), 0, 3);
    under_step = under_step + steplen;
end

figure;
imagesc(sumtvi);
colormap gray;
axis('square');

figure;
tvisum_mask = imbinarize(sumtvi, "adaptive", "Sensitivity", 0.001);
imagesc(tvisum_mask);
colormap gray;
axis('square');

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

%% Drawing script
% mask method
tviffilt=zeros(size(tvif));
tvimasks=zeros(size(tvif));
thresh=0.035;
noframes=30;
for i=1:noframes%length(frames)
	[tviffilt(:,:,i),tvimasks(:,:,i)]=maskFilter(tvif(:,:,i),thresh);
end

draw_pic(tvif(:,:,1:noframes), tvimasks, tviffilt);

%%
% medfilt2
tviffilt=zeros(size(tvif));
%thresh=0.035;
for i=1:length(frames)
	tviffilt(:,:,i)=medfilt2(tvif(:,:,i),[5,5]);
end

draw_pic(tvif, tviffilt);
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

function line_filtered = full_line_filter(tvi)
    xmax = size(tvi,1);
    ymax = size(tvi,2);
    zmax = size(tvi,3);
    line_filtered=zeros(xmax,ymax,zmax);
    for i=1:zmax
        tvitemp = mat2gray(tvi(:,:,i));
        tvi1 = imbinarize(tvitemp,'adaptive','Sensitivity',0.1);
        tvi2 = line_masker(tvi1, true);
        tvi1 = line_masker(tvi1, false);
        tvi2 = line_masker2(tvi2, tvi1);
        line_filtered(:,:,i) = regionfill(tvi(:,:,i),tvi2);
    end
end

function masked = line_masker(mask, begin)
    fillLength = 20;
    fillRatio = 0.9;
    perc = 0.05;
    masked = mask;
    xmax = size(masked,1);
    u = xmax*(1-perc);
    l = xmax*perc;

    [xi,yi] = find(masked);

    len = length(xi);
    halflen = round(len/2);
    for i=1:halflen
        xpos = xi(i);
        ypos = yi(i);
        

        if begin && ((xpos<u)&&(xpos>l))
            masked(xpos,ypos) = 0;
            continue;
        end

        clampedRange = max(1,xpos-fillLength):min(xmax,xpos+fillLength);
        clampedRangeValues = masked(clampedRange,ypos);
        
        count = sum(clampedRangeValues);
        indices = find(clampedRangeValues);
        firstIndex = min(indices);
        lastIndex = max(indices);
        
        %fprintf("count=%d, firstIndex=%d, lastIndex=%d", count, firstIndex, lastIndex);
    
        clampedRange2 = max(1,xpos-fillLength+firstIndex):min(xmax,xpos+lastIndex);
        if count / fillLength > fillRatio
            masked(clampedRange2,ypos) = 1;
        else
            %masked(clampedRange,ypos) = 0;
        end
    end
    
    for i=flip(halflen:len)
        xpos = 0;
        ypos = 0;
        try
            xpos = xi(i);
            ypos = yi(i);
        catch
            fprintf('i=%d', i);
        end

        if begin && ((xpos<u)&&(xpos>l))
            masked(xpos,ypos) = 0;
            continue;
        end

        clampedRange = max(1,xpos-fillLength):min(xmax,xpos+fillLength);
        clampedRangeValues = masked(clampedRange,ypos);
        
        count = sum(clampedRangeValues);
        indices = find(clampedRangeValues);
        firstIndex = min(indices);
        lastIndex = max(indices);
        
        %fprintf("count=%d, firstIndex=%d, lastIndex=%d", count, firstIndex, lastIndex);
    
        clampedRange2 = max(1,xpos-fillLength+firstIndex):min(xmax,xpos+lastIndex);
        if count / fillLength > fillRatio
            masked(clampedRange2,ypos) = 1;
        else
            %masked(clampedRange,ypos) = 0;
        end
    end

end

function masked = line_masker2(mask1, mask2)
    xmax = size(mask1,1);
    ymax = size(mask1,2);
    masked = zeros(xmax, ymax);

    for i=1:ymax
        if mask1(1,i) < 0.5
            continue;
        end
        
        for j=1:xmax
           if mask2(j,i) > 0.5
                masked(j,i) = 1;
           else
               break;
           end
        end

        if mask1(xmax,i) < 0.5
            continue;
        end

        for j=flip(1:xmax)
           if mask2(j,i) > 0.5
                masked(j,i) = 1;
           else
               break;
           end
        end
    end
end

function draw_pic(mat1, mat2, mat3, mat4, delay)
    lowi = 0.0;
    highi = 0.1;

    figure(1);
    Q = size(mat1,3);

    W1 = mat1(:,:,1);

    W2 = mat2(:,:,1);

	W3 = mat3(:,:,1);
	
	if mat4
		W4 = mat4(:,:,1);
		W4 = mat2gray(W4);
	end

	% Maximize window.
	g = gcf;
	g.WindowState = 'maximized';
	drawnow;

    subplot(2,2,1);
    img1 = imagesc(W1, [lowi, highi]);

    subplot(2,2,2);
    img2 = imagesc(W2, [lowi, highi]);

	subplot(2,2,3);
    img3 = imagesc(W3, [lowi, highi]);
	colorbar;

	if mat4
		subplot(2,2,4);
		img4 = imagesc(W4, [lowi, highi]);
	end

    pause(2);
    for K = 2:Q
        W1 = mat1(:,:,K);
        W1 = imadjust(W1, [0,1],[]);

        W2 = mat2(:,:,K);
        W2 = imadjust(W2, [0,1],[]);
		
		W3 = mat3(:,:,K);
        W3 = imadjust(W3, [0,1],[]);

		if mat4
			W4 = mat4(:,:,K);
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
    img1 = imshow(W1, [lowi, highi]);
    axis('square')

    subplot(1,2,2);
    img2 = imshow(W2, [lowi, highi]);
    axis('square')
    
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
