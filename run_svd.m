function svd_dat = run_svd(svd_dat,upper,lower)
    shape = size(svd_dat);
    svd_dat = reshape(svd_dat, shape(1)*shape(2), shape(3));
    
    [U,S,V] = svd(svd_dat, 'econ');
    
    S((lower+1):end, (lower+1):end) = 0;
    
    S(1:(upper-1),1:(upper-1)) = 0;

    svd_dat = single(reshape(U*S*V', shape(1), shape(2), shape(3)));
end