function stdratio = roi_metric(data, roi)
% Viktors contrast computing script based on data and roi independent of dataformat.

    stdmat = std(data,[],3);
    roistd = mean(stdmat(roi),"all");
    outsidestd = mean(stdmat(~roi), "all");
    stdratio = roistd / outsidestd;
end