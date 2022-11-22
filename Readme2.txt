% How to read this large file (in pieces)
fn_STR = '181023_1311.mat';

h = matfile(fn_STR);
% Get the "sound data" or "Radiofrequency data". 
% The columns, i.e. second dimension is the channels (x-direction) and y-diretion is the sound signa which is also depth dimension y
RF_MAT = h.rf_data_set(:,:,1:4000); % Read first 4000 frames. NB: Frames are samples at 2000Hz, i.e. 0.5ms between each frame. so 2000 frames is 1 second sequence

% How to reconstruct an image from the RF_MAT data
frame_select = 10;
B_MAT = sqrt(abs(hilbert(squeeze(RF_MAT(:,:,frame_select))))); % The sqrt is just for "compression" or the intensity data


% You can read the corresponding TVI data using the same approach as above.
