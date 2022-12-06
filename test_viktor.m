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
Snew(1:4,1:4) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%%
zmax = size(Bmodes_new,3);
steplen = round(zmax/5);
under_step = 1;
sumbmode = zeros(size(Bmodes_new,1),size(Bmodes_new,2));
for i=steplen:steplen:zmax
    sumbmode = sumbmode + std(Bmodes_new(:,:,under_step:i), 0, 3);
    under_step = under_step + steplen;
end

%%
figure;
imagesc(sumbmode);
colormap gray;
axis('square');

figure;
bmodesum_mask = imbinarize(sumbmode, "adaptive");
imagesc(bmodesum_mask);
colormap gray;
axis('square');

%%

%draw_pic(tvif, Bmodes_new);
%draw_pic2(tvi_line_filtered, Bmodes_new);
draw_pic2(tvif, tvi_line_filtered);
%draw_pic(tvi_line_filtered);

%%
tviff = medfilt3(tvi_line_filtered);

%%

draw_pic2(tvi_line_filtered, tviff);
%%

after_lf_masks = logical(zeros(size(tvi_line_filtered)));
for i=1:size(tvi_line_filtered,3)
    after_lf_masks(:,:,i) = imbinarize(tvi_line_filtered(:,:,i),"adaptive","ForegroundPolarity","dark");
end

%%
draw_pic2(tviff, after_lf_masks);

%%
tviffconv = zeros(size(tviff));
for i=1:size(tviff,3)
    tviffconv(:,:,i) = conv2(tviff(:,:,i),ones(4,2)/8,'same');
end
%%
draw_pic2(tviff, tviffconv);

%%
stdpic = std(tviffconv,0,3);

%%
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
draw_pic2(tviffconv, tlf_pf);

%%
tlen = size(tvi_line_filtered, 3);
zmax = size(tvi_line_filtered, 3);
xstart = 1200;
xstop = 1300;
ystart = 80;
ystop = 128;

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

