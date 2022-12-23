function draw_pic2(mat1, mat2, pause_time, perc)
    if nargin < 3
        pause_time = 0.02;
        perc = [1, 99];
    elseif nargin < 4
        perc = [1, 99];
    end

    ptiles1 = prctile(mat1,perc,[1,2,3]);
    ptiles2 = prctile(mat2,perc,[1,2,3]);

    figure(1);
    Q = size(mat1,3);
    W1 = mat1(:,:,1);

    W2 = mat2(:,:,1);

    subplot(1,2,1);
    img1 = imagesc(W1);
    clim(ptiles1);

    subplot(1,2,2);
    img2 = imagesc(W2);
    clim(ptiles2);
    
    for K = 2:Q
        W1 = mat1(:,:,K);
        W2 = mat2(:,:,K);
        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
        %drawnow limitrate;
        drawnow();
        pause(pause_time);
    end
end