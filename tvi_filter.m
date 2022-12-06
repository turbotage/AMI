function tvi_dat = tvi_filter(tvi_dat, do_fft_peak_removal, ...
    do_butter, do_hardline, do_whole_vol, do_slice_conv)
    
    % Frequency peak removal
    if do_fft_peak_removal
        zmax = size(tvi_dat,3);
        ftvi = fft(tvi_dat,[],3);
        %rmid = 33;
        %ftvi(:,:,(rmid-1):(rmid+1)) = 0; ftvi(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
        %ftvi(:,:,rmid) = 0; ftvi(:,:,zmax-rmid+2) = 0;
        %rmid = 17;
        %ftvi(:,:,(rmid-1):(rmid+1)) = 0; ftvi(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
        %ftvi(:,:,rmid) = 0; ftvi(:,:,zmax-rmid+2) = 0;
        rmid = 49;
        %ftvi(:,:,(rmid-1):(rmid+1)) = 0; ftvi(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
        ftvi(:,:,rmid) = 0; ftvi(:,:,zmax-rmid+2) = 0;

        tvi_dat = ifft(ftvi,[],3);
    end

    % Lowpass in temporal
    if do_butter
        [b,a] = butter(1, 200/250,'low');
        tvi_dat = filter(b,a,tvi_dat,[],3);
    end

    % Hard line removal
    if do_hardline
        tvi_dat = full_line_filter(tvi_dat);
    end

    % General whole volume filter
    if do_whole_vol
        tvi_dat = medfilt3(tvi_dat);
    end

    % Per temporal slice convolution
    if do_slice_conv
        for i=1:size(tvi_dat,3)
            tvi_dat(:,:,i) = conv2(tvi_dat(:,:,i),ones(3,2)/6,'same');
        end
    end

end