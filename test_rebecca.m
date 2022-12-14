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
frames=1:250;
RF_MAT= single(squeeze(RF_MAT(:,:,frames)));
Bmodes = single(sqrt(abs(hilbert(squeeze(RF_MAT(:,:,frames))))));
shape = size(Bmodes);
RF_MAT_f=reshape(RF_MAT, shape(1)*shape(2), shape(3));
%% svd
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(25:end, 25:end) = 0;
Snew(1:20,1:20) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%% H-scan
% generate kernel of hermitian polynomials of ord 2 and 8
% with matrices
b=6;
a=20;
fs=1000;
t = -5:1/fs:5;
GH2=exp(-(t/b).^2).^2.*(4*t.^2-2).*exp(-t.^2);
GH8=exp(-(t/a).^2).^2.*((2*t).^8-3584*t.^6+13440*(t.^4-t.^2)+1680).*exp(-t.^2);
E2 = 1*3*sqrt(pi/2);
E8 = 1*3*5*7*9*11*13*15*sqrt(pi/2);

figure
sgtitle('Gaussian weighted Hermite polynomials')
subplot(2,1,1)
plot(t,GH2,'.b')
title('GH_2')
ylim([min(GH2,[],'all')*1.25 max(GH2,[],'all')*1.25]);
subplot(2,1,2)
plot(t,GH8,'.r')
title('GH_8')
ylim([min(GH8,[],'all')*1.25 max(GH8,[],'all')*1.25]);

%% fft (envelopes) of original RF signal and filters
clear B0 B2 B8
frame=90;

[pxx,f]=pwelch(RF_MAT_f(:,frame));
[pxx2,f]=pwelch(GH2);
[pxx8,f]=pwelch(GH8);

figure
plot(f,pxx,'color',[0 .5 0],'linestyle','-');
hold on; plot(f,pxx8,'-b');
hold on; plot(f,pxx2,'-r');
legend('Original RF-signal','High pass (GH8)','Low pass (GH2)');
xlabel('Frequency f [Hz]'); ylabel('Power per frequency P [dB/Hz]');
sgtitle(['PSD for frame ',num2str(frame)])
legend('Original RF-signal','High pass (GH8)','Low pass (GH2)');
hold off
hold off


%% RGB code signal and convolute
clear RF_MAT2 RF_MAT8 
clear Bmodes2 Bmodes8

% convolution
noframes=5;
for j=1:noframes
	RF_MAT2(:,:,j)=conv2(GH2,RF_MAT(:,:,j))./sqrt(E2);
	RF_MAT8(:,:,j)=conv2(GH8,RF_MAT(:,:,j))./sqrt(E8);
	j
end

% plot H-scan filtering results
Bmodes2=sqrt(abs(hilbert(RF_MAT2)));
Bmodes8=sqrt(abs(hilbert(RF_MAT8)));

draw_pic(Bmodes(:,:,1:noframes), Bmodes2, Bmodes8, 0.03);

%% combination with colorcoding in freq space to RGB and then combined into 
% a new image RGB colormap R for low freqs, B for high, G for original sig
% envelope.

% RGB channeling intensites G=signal, B=GH8, R=GH2
% mymap = [	1 0 0 ;
%			0 1 0 ;
%			0 0 1];
% ind2rgb(BmodesTot,mymap)

%% load TVI (Validation) data
fn_STR = '181023_1311tvi_rs.mat';
h = matfile(fn_STR);
tvif_d = single(h.tvi_downsampled(:,:,1:250));
tvif_d = sqrt(abs(hilbert(squeeze(tvif_d(:,:,:)))));

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

[b,a] = butter(1, 200/250,'low');
tvif = filter(b,a,tvif_d,[],3);
Bmodes_new = filter(b,a,Bmodes_new,[],3);

for j=1:noframes
	tviff(:,:,j) = medfilt2(tvi_line_filtered(:,:,j));
	Bmodes_newff(:,:,j) = medfilt2(Bmodes_new(:,:,j));
end

%%

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
[b,a] = butter(4,[5 25]./(Fsamp/2),'bandpass');
y_filt = filtfilt(b,a,y);

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

function draw_pic(mat1, mat2, mat3, delay)
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

	% Maximize window.
	g = gcf;
	g.WindowState = 'maximized';
	drawnow;

    subplot(1,3,1);
    img1 = imshow(W1, [lowi, highi]);
    axis('square')

    subplot(1,3,2);
    img2 = imshow(W2, [lowi, highi]);
    axis('square')

	subplot(1,3,3);
    img3 = imshow(W3, [lowi, highi]);
    axis('square')

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
        caxis([0,1])
        pause(delay);
        disp(K);
    end
end

function draw_pic2(mat1, mat2)
    lowi = 0.0;
    highi = 0.1;

    figure(1);
    Q = size(mat1,3);
    W1 = mat1(:,:,1);
    W1 = mat2gray(W1);
    %W1 = imadjust(W1, [0,1],[]);

    W2 = mat2(:,:,2);
    W2 = mat2gray(W2);
    %W2 = imadjust(W2, [0,1],[]);


    subplot(1,2,1);
    img1 = imshow(W1, [lowi, highi]);
    axis('square')

    subplot(1,2,2);
    img2 = imshow(W2, [lowi, highi]);
    axis('square')
    
    pause(2);
    for K = 2:Q
        W1 = mat1(:,:,K);
        W1 = mat2gray(W1);
        W1 = imadjust(W1, [0,1],[]);
        W2 = mat2(:,:,K);
        W2 = mat2gray(W2);
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
