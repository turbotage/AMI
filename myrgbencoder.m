function encoded_channel = myrgbencoder(data)
% Rebeccas script which encodes each frame of the data's intensities to a 
% value between 0 and 1 saving it and returning it as one rgb channel matrix. 
% It can handle 2D matrices both with and without a temporal dimension.

    if size(data,3)
        encoded_channel = single(zeros([size(data,1) size(data,2) 1 size(data,3)]));
        for j=1:size(data,3)
            encoded_channel(:,:,1,j) =  (data(:,:,j)-min(data(:,:,j),[],'all')) / ...
                                        (max(data(:,:,j),[],'all')-min(data(:,:,j), ...
                                        [],'all'));
        end
        
    else
        encoded_channel = single(zeros(size(data)));
        encoded_channel = (data(:,:)-min(data(:,:),[],'all')) / ...
                          (max(data(:,:),[],'all')-min(data(:,:),[],'all'));
    end
end