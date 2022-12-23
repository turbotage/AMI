function stdratio = roi_metric(data, roi)
    stdmat = std(data,[],3);
    roistd = mean(stdmat(roi),"all");
    outsidestd = mean(stdmat(~roi), "all");
    stdratio = roistd / outsidestd;
end