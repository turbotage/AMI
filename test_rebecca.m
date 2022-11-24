%  AMI Project 2022
%% load dataset 1 ... only image contraction w ~ Hz full activation vs full rest

%% load dataset 2 stimulated contraction w ~1.6Hz 1 muscle complex
fn_STR = '181023_1311.mat';
h = matfile(fn_STR);
Fsamp=2000;
rfmatf = single(h.rf_data_set(:,:,1:500)); % need full set but take every 
% other or every fourth, downsampling to 1000 or 500 Hz

%% load dataset 3 voluntary contraction w ~5-25Hz 2-4 muscle complexes
Bmodes = sqrt(abs(hilbert(squeeze(rfmatf(:,:,:)))));
clear rfmatf
shape = size(Bmodes);
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
clear Bmodes

%% load TVI data
fn_STR = '181023_1311tvi.mat';
h = matfile(fn_STR);
tvif = single(h.TVI_MAT(:,:,1:500));
tvif = sqrt(abs(hilbert(squeeze(tvif(:,:,:)))));

%% SVD
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(2:end, 25:end) = 0;
Snew(1:4,1:4) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%% SVD vs TVI (Validation)
% Prefiltering
% median3 and butter order 4 on both TVI and SVD, keeping the freqs in
% range 5-25 Hz

% tripple median filtering
Bmodes_fin=medfilt3(Bmodes_new);
tvi_fin=medfilt3(tvif);

% spatial filtering
C = conv2(A, B);

% freq filtering
% idea: bandpass on every pixelintensity over activationtime 5-15(30) Hz
% for dataset 2 but around 1.5 Hz for dataset 3
[b,a] = butter(4,[5 25]./(Fsamp/2),'bandpass');
y_filt = filtfilt(b,a,y);

%% H-scan
% generate kernel of hermitian polynomials of ord 2 and 8 w gaussian 
% windowing

% convolution

% combination with colorcoding in freq space to RGB and then combined into 
% a new image RGB colormap R for low freqs, B for high, G for original sig
% envelope.


%% H-scan vs TVI (Validation)
% Prefiltering
% median3 and butter order 4 on both TVI and SVD, keeping the freqs in
% range 5-25 Hz

% tripple median filtering
Bmodes_fin=medfilt3(Bmodes_new);
tvi_fin=medfilt3(tvif);

% spatial filtering
C = conv2(A, B);

% freq filtering
% idea: bandpass on every pixelintensity over activationtime 5-15(30) Hz
% for dataset 2 but around 1.5 Hz for dataset 3
[b,a] = butter(4,[5 25]./(Fsamp/2),'bandpass');
y_filt = filtfilt(b,a,y);

%% Drawing script

draw_pic(tvif, Bmodes_new);

function draw_pic(mat1, mat2)
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