function draw_pic2(mat1, mat2)
    figure(1);
    Q = size(mat1,3);
    W1 = mat1(:,:,1);
    W1 = mat2gray(W1);

    W2 = mat2(:,:,1);
    W2 = mat2gray(W2);

    subplot(1,2,1);
    img1 = imshow(W1);
    axis('square')

    subplot(1,2,2);
    img2 = imshow(W2);
    axis('square')
    
    pause(2);
    for K = 2:Q
        W1 = mat1(:,:,K);
        W1 = mat2gray(W1);
        W2 = mat2(:,:,K);
        W2 = mat2gray(W2);
        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
        %drawnow limitrate;
        drawnow();
        pause(0.02);
    end
end