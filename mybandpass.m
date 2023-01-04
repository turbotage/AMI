function bandpassed = mybandpass(data, order, cutoffs, Fs)
    [b,a] = butter(order, cutoffs / (Fs/2), "bandpass");
%     freqz(b,a);
%     keyboard
    bandpassed = single(zeros(size(data)));
    for i=1:size(data,1)
        for j=1:size(data,2)
            bandpassed(i,j,:) = single(filtfilt(b,a,double(squeeze(data(i,j,:)))));
        end
    end
end