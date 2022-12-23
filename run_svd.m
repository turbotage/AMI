function svd_dat = run_svd(svd_dat,upper,lower)
    shape = size(svd_dat);
    svd_dat = reshape(svd_dat, shape(1)*shape(2), shape(3));
    
    [U,S,V] = svd(svd_dat, 'econ');
    
    if lower <= size(svd_dat,2)
        S(lower:end, lower:end) = 0;
    end
    
    if upper > 0
        S(1:upper,1:upper) = 0;
    end

    svd_dat = single(reshape(U*S*V', shape(1), shape(2), shape(3)));
end