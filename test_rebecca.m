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
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));

%% svd
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(25:end, 25:end) = 0;
Snew(1:20,1:20) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%% H-scan
% generate kernel of hermitian polynomials of ord 2 and 8
clear RF_MAT2 Bmodes2 RF_MAT8 Bmodes8

RF_MAT=reshape(RF_MAT, shape(1)*shape(2), shape(3));

noframes=250;

% HOW WE SHOULD WE GET THE GAUSSIAN WEIGHTS? (SEEMS WIERD AS IT
% DOESN'T CHANGE THE RESULTING POLYNOMIAL MUCH.
% HOW DO WE GET THE FILTERS A AND B (CENTER FREQ) SETTINGS? AND WHAT IS
% TAU?

% Windowing?
syms x y
H2=hermiteH(2,x);
H8=hermiteH(8,x);

H2=sym2poly(H2);
H8=sym2poly(H8);

% convolution
E2 = 1*3*sqrt(pi/2);
E8 = 1*3*5*7*9*11*13*15*sqrt(pi/2);

[H2,~] = smoothdata(H2,'gaussian');
[H8,~] = smoothdata(H8,'gaussian');

RF_MAT2=conv2(H2,RF_MAT)./sqrt(E2);
RF_MAT8=conv2(H8,RF_MAT)./sqrt(E8);

%% fft (envelopes) of original RF signal and filters
clear B0 B2 B8
frame=90;

% HOW DO WE GET THE FILTERS ENVELOPES?
% muH2=mean(B2);
% stdH2=std(B2);
% x = 1:2e5;
% GH2 = exp(- 0.5 * ((x - muH2) / stdH2) .^ 2) / (stdH2 * sqrt(2 * pi));
% muH8=mean(B8);
% stdH8=std(B8);
% GH8 = exp(- 0.5 * ((x - muH8) / stdH8) .^ 2) / (stdH8 * sqrt(2 * pi));

% WHAT IS THE SAMPLING FREQ FOR RF_MAT? (500 imgs of a muscl contraction.)
fs=1000;
figure
plot(fft(RF_MAT(:,frame)));hold on;
plot(fft(RF_MAT2(:,frame)));hold on;
plot(fft(RF_MAT8(:,frame)))
legend('Original RF-signal','High pass (GH8)','Low pass (GH2)');
xlabel('Frequency f [Hz]');ylabel('Power per frequency P [dB/Hz]');
hold off

[B0]=envelope(RF_MAT(:,frame),2000,'peak');
[B2]=envelope(RF_MAT2(:,frame),2000,'peak');
[B8]=envelope(RF_MAT8(:,frame),2000,'peak');
% B0=sqrt(abs(B0));
% B2=sqrt(abs(B2));
% B8=sqrt(abs(B8));
% Filter overview plot
figure
plot(B0,'color',[0 .5 0],'linestyle','-'); hold on;
plot(B2,'-r');hold on;
plot(B8,'-b');
sgtitle(['Envelopes for frame ',num2str(frame)])
legend('Original RF-signal','High pass (GH8)','Low pass (GH2)');
xlabel('Time t [s]');ylabel('Amplitude A [1]');
hold off

%% RGB code signal

%% plot H-scan results
Bmodes2=reshape(RF_MAT2, shape(1), shape(2), noframes+2);
Bmodes2=sqrt(abs(hilbert(Bmodes2)));

Bmodes8=reshape(RF_MAT8, shape(1), shape(2), noframes+8);
Bmodes8=sqrt(abs(hilbert(Bmodes8)));

draw_pic(Bmodes(:,:,1:noframes), Bmodes2, Bmodes8);

% combination with colorcoding in freq space to RGB and then combined into 
% a new image RGB colormap R for low freqs, B for high, G for original sig
% envelope.

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
function [x, w] = GaussHermite_2(n)
% This function determines the abscisas (x) and weights (w) for the
% Gauss-Hermite quadrature of order n>1, on the interval [-INF, +INF].
    % This function is valid for any degree n>=2, as the companion matrix
    % (of the n'th degree Hermite polynomial) is constructed as a
    % symmetrical matrix, guaranteeing that all the eigenvalues (roots)
    % will be real.
    
	% © Geert Van Damme
	% geert@vandamme-iliano.be
	% February 21, 2010    

	% Building the companion matrix CM
    	% CM is such that det(xI-CM)=L_n(x), with L_n the Hermite polynomial
    	% under consideration. Moreover, CM will be constructed in such a way
    	% that it is symmetrical.
		i   = 1:n-1;
		a   = sqrt(i/2);
		CM  = diag(a,1) + diag(a,-1);

	% Determining the abscissas (x) and weights (w)
    	% - since det(xI-CM)=L_n(x), the abscissas are the roots of the
    	%   characteristic polynomial, i.d. the eigenvalues of CM;
    	% - the weights can be derived from the corresponding eigenvectors.
		[V L]   = eig(CM);
		[x ind] = sort(diag(L));
		V       = V(:,ind)';
		w       = sqrt(pi) * V(:,1).^2;
end

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

function draw_pic(mat1, mat2, mat3)
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
        pause(0.03);
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
