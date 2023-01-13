function datamat = svd_sweep(data,roi,sstart,inters)
    datamat = zeros([length(sstart), length(inters)]);
    tlen = size(data,3);

    for i=1:length(sstart)
        for j=1:length(inters)
            start = round(tlen*sstart(i));
            start = max([start, 1]);
            inter = inters(j);
            svd_dat = run_svd(data,start,start+inter);

            ratio = roi_metric(svd_dat, roi);

            datamat(i,j) = ratio;

            fprintf("%d of %d starts, %d of %d inters\n", i,length(sstart),j,length(inters));
        end
    end
end