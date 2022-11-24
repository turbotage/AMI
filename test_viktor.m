%%
fn_STR = '181023_1311.mat';
ndata = 500
h = matfile(fn_STR);
rfmat = h.rf_data_set(:,:,1:4000);
rfmat_downsampled = resample(rfmat,1, 4,'Dimension', 3);

%%


%%
rfmat_dsf = single(load('181023_1311_rs.mat').rfmat_downsampled);

%%
%load DATA
Bmodes = sqrt(abs(hilbert(squeeze(rfmat_dsf(:,:,:)))));
shape = size(Bmodes)
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
clear Bmodes

%%
% load TVI
fn_STR = '181023_1311tvi.mat';
h = matfile(fn_STR);
tvif = single(h.TVI_MAT(:,:,1:500));
tvif = sqrt(abs(hilbert(squeeze(tvif(:,:,:)))));

%%
[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(25:end, 25:end) = 0;
Snew(1:20,1:20) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%figure(1);
%animate_stuff(Bmodes);
%figure(2);
%animate_stuff(Bmodes_new);
%%

draw_pic(tvif, Bmodes_new);

function draw_pic(mat1, mat2)
    lowi = 0.0
    highi = 0.1

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
