new_viktor.m

Runs a parametric sweep for singular value intervals and produces an std map for tvi and svd
for a chosen singular value interval.

Choose dataset via the strings in the load sections, 
choose bandpass appropriatly for the specific dataset.


tvi_contrasts.m

Calculates the contrast for the TVI method for a specific dataset using roi_metric

svd_contrasts.m

Calculates the contrast for a specific dataset with the SVD method using roi_metric,
specific singular value interval too.

bp_contrasts.m

Same as svd_contrasts.m and tvi_contrasts.m but for just a bandpass


svd_data2_plots.m

Plots the contrast surface plot for a parametric sweep. 
Ratios matrix should be stored from parametric sweep

resample_data.m

Takes in ultrasound datasets, resamples them to 1/4 of the original sampling rate and stores
them as single. Also, only the upper 1000 voxels are stored. This often decreases the storage
by a factor 10-16.