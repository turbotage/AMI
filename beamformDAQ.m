function out = beamformDAQ(preBeamformed, fulAprSz)
% The following script take the pre-beamform data and creates beamfromed
% image by applying parallel beamformers to the same data. 
%
% input:
%   preBeamformed:  the first dimension should be scan lines and the second
%                   dimension should be the channel samples
%   fulAprSz:       size of the aperture beamformer
%
% output:
%   out:            beamformed data
%
% Example:
%   bRF = beamformDAQ(RF, 32);
%
% Author: Reza Zahiri Azar
% Copyright Ultrasonix Medical Corporation - April 2010

% The first dimension has to be lines and second dimension has to be samples
if (size(preBeamformed,1) > size(preBeamformed,2))
    preBeamformed = preBeamformed';
    transposed = true;
else
    transposed = false;
end

nl = size(preBeamformed,1);   % number of lines
ns = size(preBeamformed,2);   % number of samples

hlfAprSz = floor(fulAprSz/2);
fulAprSz = hlfAprSz*2 + 1;

win = hanning(fulAprSz)';   % windowing used for apodization
chnls = zeros(1,fulAprSz);  % will be filled with proper signal from each channel

channelSpacing = 0.3;    % spacing between lines/channels in the lateral direction in milimeter
sampleSpacing = 0.01925; % spacing between samples in the axial direction in milimeter

offset = 10; %10;             % distance between first sample and transducer element in milimeter
nSampleOffset = round(offset/ sampleSpacing);    % offset in samples: 1st sample is actually 1+nSampleOffset sample

x = -hlfAprSz: hlfAprSz;    % aperture indices

for i = 1:nl    % for each line/channel
%     disp(['scanline #', num2str(i)]);
    
    for j=1:ns  % find the value for each sample
        
        % calc delays based on sample depth and receive aperture
        distance2probe =  (j + nSampleOffset) * sampleSpacing ; 
        rad = sqrt( (hlfAprSz*channelSpacing)^2 +  distance2probe^2 );
        scale = (rad - distance2probe);
        delays = + x.^2 / (hlfAprSz)^2 * scale/ sampleSpacing  ; 

        cntr = i;       % center of apreture
        apr = cntr + x; % full aperture indices
        
        chnls = zeros(1,fulAprSz);  % initialize channel values
        cnt = 0;            % count the number of channel that we use
                            % this parameter can be used for averaging (not used for now)

        % find the corresponding value from each channel
        for k = 1:fulAprSz  
            
            chlIndx = apr(k);
            if chlIndx<1, continue, end;
            if chlIndx>nl, continue, end;
            cnt = cnt + 1;

            chlSmpl = round(j+delays(k));

            if chlSmpl<1, continue, end;
            if chlSmpl>ns, continue, end;
            
            chnls(k) = preBeamformed(chlIndx, chlSmpl); 
            
        end;
        % apodization : ideally has to be a function of depth
        chnls = win .* chnls;
        
        % beamfroming
        out(i,j) = sum( chnls ) ;%/ cnt;
        
    end;
end
if (transposed),    out = out'; end

