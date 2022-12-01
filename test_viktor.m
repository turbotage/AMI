% load DATA
rfmat_dsf = single(load('181023_1311_rs.mat').rfmat_downsampled);

%%
%load DATA
Bmodes = sqrt(abs(hilbert(squeeze(rfmat_dsf(:,:,:)))));
shape = size(Bmodes)
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
clear Bmodes

%%
% load TVI
fn_STR = '181023_1311tvi_rs.mat';
h = matfile(fn_STR);
tvif_d = single(h.tvi_downsampled(:,:,1:1000));

%%
tvif_d = sqrt(abs(hilbert(squeeze(tvif_d(:,:,:)))));

%%

[b,a] = butter(1, 200/250,'low');
tvif = filter(b,a,tvif_d,[],3);

%%
lowi = 0;
highi = 1;
frame = 51;

tvi1 = mat2gray(tvif(:,:,frame));
tvi1a = imbinarize(tvi1,'adaptive','Sensitivity',0.001);
tvi1b = line_masker(tvi1a, true);
tvi1a = line_masker(tvi1a, false);
tvi1b = line_masker2(tvi1b, tvi1a);

tvi2 = regionfill(tvi1,tvi1b);

%tvi1e = and(tvi1a,tvi1b);
figure;
subplot(2,1,1);
imagesc(tvi1a);
axis('square');
caxis([0, 1]);
title('tvi1a')

subplot(2,1,2);
imagesc(tvi1b);
axis('square');
caxis([0, 1]);
title('tvi1b')

figure;
subplot(2,1,1);
imagesc(tvi1);
axis('square');
caxis([0, 1]);
title('tvi1')

subplot(2,1,2);
imagesc(tvi2);
axis('square');
caxis([0, 1]);
title('tvi2');

%%
tvi_line_filtered = full_line_filter(tvif);

%%
figure;
imagesc(tvi_line_filtered(:,:,1));
axis('square');
caxis([0,1])

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

%draw_pic(tvif, Bmodes_new);
draw_pic2(tvi_line_filtered, Bmodes_new);


