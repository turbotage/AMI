function draw_pic2(mat1, mat2)
    lowi = 0.0;
    highi = 0.1;

    figure(1);
    Q = size(mat1,3);
    W1 = mat1(:,:,1);
    W1 = mat2gray(W1);
    %W1 = imadjust(W1, [0,1],[]);

    W2 = mat2(:,:,2);
    W2 = mat2gray(W2);
    %W2 = imadjust(W2, [0,1],[]);


    subplot(1,2,1);
    img1 = imshow(W1, [lowi, highi]);
    axis('square')

    subplot(1,2,2);
    img2 = imshow(W2, [lowi, highi]);
    axis('square')
    
    pause(2);
    for K = 2:Q
        W1 = mat1(:,:,K);
        W1 = mat2gray(W1);
        W1 = imadjust(W1, [0,1],[]);
        W2 = mat2(:,:,K);
        W2 = mat2gray(W2);
        W2 = imadjust(W2, [0,1],[]);
        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
        %drawnow limitrate;
        drawnow();
        caxis([0,1])
        pause(0.02);
        disp(K);
    end
end