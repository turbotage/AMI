%% load TVI
tvi_dat = matfile('181023_1555tvi_rsf.mat').tvi_dat;
disp('Ran load TVI')
% bandpass [3,20] for 1311 and 1555, [0.5,10] for 0913
bp_interval = [3,20];
%% load ROI
roi = matfile('181023_1555_ROI.mat').ROI;
disp('Ran load roi')
%% TVI Filter
tvi_dat_f = mybandpass(tvi_dat,4,bp_interval,500);
tvi_dat_f = mymedfilt(tvi_dat_f, [10,3]);
%tvi_dat_f = full_line_filter(tvi_dat_f);
disp('TVI Filter')
%%
contrast = roi_metric(tvi_dat_f,roi);
