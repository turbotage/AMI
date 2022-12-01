% load DATA
rfmat_dsf = single(load('181023_1311_rs.mat').rfmat_downsampled);

%%
%load DATA
Bmodes = sqrt(abs(hilbert(squeeze(rfmat_dsf(:,:,:)))));
shape = size(Bmodes)
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));
clear Bmodes
clear rfmat_dsf

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
tvi_line_filtered = full_line_filter(tvif);

%%
figure;
imagesc(tvi_line_filtered(:,:,51));
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
%draw_pic2(tvi_line_filtered, Bmodes_new);
draw_pic2(tvif, tvi_line_filtered);
%draw_pic(tvi_line_filtered);

%%
tviff = medfilt3(tvi_line_filtered);

%%
draw_pic2(tvif, tviff);

%%
stdpic = std(tviff,0,3);

%%
stdpic = conv2(stdpic,ones(5,5)/25, 'same');

%%
figure;
imagesc(stdpic);
colormap gray;
axis('square');



%%
tlen = size(tvi_line_filtered, 3);
xstart = 1200;
xstop = 1300;
ystart = 80;
ystop = 128;

avg_intense = reshape(mean(mean( ...
    tvi_line_filtered(xstart:xstop,ystart:ystop,:),1),2),[tlen,1]);
figure;
plot(avg_intense);
freq_int = fft(avg_intense);
P2 = abs(freq_int/tlen);
P1 = P2(1:tlen/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Fs = 500;

figure;
f = Fs*(0:(tlen/2))/tlen;
plot(f,P1) 
