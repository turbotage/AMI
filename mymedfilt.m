function filtered = mymedfilt(data,neigh)
    
    if nargin < 2
        neigh = [3,3];
    end

    if size(data,4)
        filtered = single(zeros([size(data,1) size(data,2) 1 size(data,4)]));
        if any(data(:,:,1,:),'all')
            for i=1:size(data,4)
                filtered(:,:,1,i) = single(medfilt2(data(:,:,1,i),neigh));
            end

        elseif any(data(:,:,3,:),'all')
            for i=1:size(data,4)
                filtered(:,:,3,i) = single(medfilt2(data(:,:,3,i),neigh));
            end

        else
            for i=1:size(data,4)
                filtered(:,:,2,i) = single(medfilt2(data(:,:,2,i),neigh));
            end 
        end
    else
        filtered = single(zeros(size(data)));
        for i=1:size(data,3)
            filtered(:,:,i) = single(medfilt2(data(:,:,i),neigh));
        end
    end
end