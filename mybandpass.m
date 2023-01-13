function bandpassed = mybandpass(data, order, cutoffs, Fs)
% Viktors band pass filtering script using a digital iir butterworth filter of 
% specified order, cut-offs with corresponding sampling frequency for data.

    [b,a] = butter(order, cutoffs / (Fs/2), "bandpass");

%     freqz(b,a); % Previouly used to visualize the filter
%     keyboard 

    bandpassed = single(zeros(size(data)));
    for i=1:size(data,1)
        for j=1:size(data,2)
            bandpassed(i,j,:) = single(filtfilt(b,a,double(squeeze(data(i,j,:)))));
        end
    end
end