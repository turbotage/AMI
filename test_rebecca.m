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
RF_MAT=single(squeeze(RF_MAT(:,:,frames)));
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
% generate kernel of hermitian polynomials of ord 2 and 8 w gaussian 
% windowing
syms x2 x8
H2=sym2poly(hermiteH(2, x2));
H8=sym2poly(hermiteH(8, x8));

% IS THIS HOW WE SHOULD GET THE GAUSSIAN WEIGHTS? (SEEMS WIERD AS IT
% DOESN'T CHANGE THE RESULTING POLYNOMIAL MUCH.
% HOW DO WE GET THE FILTERS A AND B (CENTER FREQ) SETTINGS? AND WHAT IS
% TAU?
[H2,window2] = smoothdata(H2,'gaussian');
[H8,window8] = smoothdata(H8,'gaussian');

clear RF_MAT2 Bmodes2 RF_MAT8 Bmodes8

RF_MAT=reshape(RF_MAT, shape(1)*shape(2), shape(3));

noframes=250;

for j=1:noframes % length(frames)
	RF_MAT(:,j) = RF_MAT(:,j)/sqrt(sum(RF_MAT(:,j)));
	RF_MAT2(:,j) = RF_MAT(:,j)/sqrt(sum(RF_MAT(:,j).^2));
	RF_MAT8(:,j) = RF_MAT(:,j)/sqrt(sum(RF_MAT(:,j).^8));
end

% convolution
RF_MAT2=conv2(H2,RF_MAT2);
RF_MAT8=conv2(H8,RF_MAT8);

% fft (envelopes) of original RF signal and filters
clear B0 B2 B8
frame=90;

figure
B0=envelope(RF_MAT(:,frame));

% HOW DO WE GET THE FILTERS ENVELOPES?
B2=envelope(polyval(H2,linspace(1,100)));
B8=envelope(polyval(H8,linspace(1,100)));

% Filter overview plot
plot(B0,'color',[0 .5 0],'linestyle','-'); hold on;
plot(B2,'-r');hold on;
plot(B8,'-b');
sgtitle(['Envelopes for frame ',num2str(frame)])
legend('Original RF-signal','High pass (GH8)','Low pass (GH2)')
hold off

%%
Bmodes2=reshape(RF_MAT2, shape(1), shape(2), noframes+2);
Bmodes2=sqrt(abs(hilbert(Bmodes2)));

Bmodes8=reshape(RF_MAT8, shape(1), shape(2), noframes+8);
Bmodes8=sqrt(abs(hilbert(Bmodes2)));

draw_pic(Bmodes(:,:,1:noframes), Bmodes2, Bmodes8);

% combination with colorcoding in freq space to RGB and then combined into 
% a new image RGB colormap R for low freqs, B for high, G for original sig
% envelope.

%% load TVI (Validation) data
fn_STR = '181023_1311tvi_rs.mat';
h = matfile(fn_STR);
tvif = single(h.tvi_downsampled(:,:,1:250));
tvif = sqrt(abs(hilbert(squeeze(tvif(:,:,:)))));

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
% idea: bandpass on every pixelintensity over activationtime 5-15(30) Hz
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
        pause(0.05);
        disp(K);
    end
end
