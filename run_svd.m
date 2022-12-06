function svd_dat = run_svd(svd_dat,upper,lower)
    shape = size(svd_dat);
    svd_dat = reshape(svd_dat, shape(1)*shape(2), shape(3));
    
    [U,S,V] = svd(svd_dat, 'econ');

    lower = floor((1-lower)*shape(3));
    S(lower:end, lower:end) = 0;
    upper = ceil(upper*shape(3));
    S(1:upper,1:upper) = 0;

    svd_dat = reshape(U*S*V', shape(1), shape(2), shape(3));
end