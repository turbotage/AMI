function draw_pic(mat1)
    lowi = 0.0;
    highi = 0.1;

    figure(1);
    Q = size(mat1,3);
    W1 = mat1(:,:,1);
    W1 = mat2gray(W1);

    img = imshow(W1);
    axis('square');
    
    pause(2);
    for K = 2:Q
        W1 = mat1(:,:,K);
        W1 = mat2gray(W1);
        %W1 = imadjust(W1, [0,1],[]);
        set(img, 'CData', W1);
        %drawnow limitrate;
        drawnow();
        %caxis([0,1])
        pause(0.02);
        %disp(K);
    end
end