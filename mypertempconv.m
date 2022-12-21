function conved = mypertempconv(data, kernel)
    conved = single(zeros(size(data)));
    for i=1:size(data,3)
        conved(:,:,i) = single(conv2(data(:,:,i), kernel,'same'));
    end
end