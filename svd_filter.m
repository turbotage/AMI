function svd_dat = svd_filter(svd_dat, do_fft_peak_removal, ...
    do_butter, do_hardline, do_whole_vol, do_slice_conv)
    
    % Frequency peak removal
    if do_fft_peak_removal
        zmax = size(svd_dat,3);
        fsvd = fft(svd_dat,[],3);
        %rmid = 33;
        %ftvi(:,:,(rmid-1):(rmid+1)) = 0; ftvi(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
        %ftvi(:,:,rmid) = 0; ftvi(:,:,zmax-rmid+2) = 0;
        %rmid = 17;
        %ftvi(:,:,(rmid-1):(rmid+1)) = 0; ftvi(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
        %ftvi(:,:,rmid) = 0; ftvi(:,:,zmax-rmid+2) = 0;
        rmid = 49;
        %ftvi(:,:,(rmid-1):(rmid+1)) = 0; ftvi(:,:,(zmax-rmid+1):(zmax-rmid+3)) = 0;
        fsvd(:,:,rmid) = 0; fsvd(:,:,zmax-rmid+2) = 0;

        svd_dat = ifft(fsvd,[],3);
    end

    % Lowpass in temporal
    if do_butter
        [b,a] = butter(1, 200/250,'low');
        svd_dat = filter(b,a,svd_dat,[],3);
    end

    % Hard line removal
    if do_hardline
        svd_dat = full_line_filter(svd_dat);
    end

    % General whole volume filter
    if do_whole_vol
        svd_dat = medfilt3(svd_dat);
    end

    % Per temporal slice convolution
    if do_slice_conv
        for i=1:size(svd_dat,3)
            svd_dat(:,:,i) = conv2(svd_dat(:,:,i),ones(3,2)/6,'same');
        end
    end

end