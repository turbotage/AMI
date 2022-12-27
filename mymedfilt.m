function filtered = mymedfilt(data,neigh)
    if nargin < 2
        neigh = [3,3];
    end
    filtered = single(zeros(size(data)));
    for i=1:size(data,3)
        filtered(:,:,i) = single(medfilt2(data(:,:,i),neigh));
    end
end