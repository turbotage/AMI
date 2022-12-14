function conved = mypertempconv(data, kernel)
    conved = zeros(size(data));
    for i=1:size(data,3)
        conved(:,:,i) = conv2(data(:,:,i), kernel,'same');
    end
end