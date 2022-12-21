function draw_pic(mat1, pause_time, perc)
    if nargin < 2
        pause_time = 0.02;
        perc = [1, 99];
    elseif nargin < 3
        perc = [1, 99];
    end

    ptiles = prctile(mat1,perc,[1,2,3]);

    figure(1);
    Q = size(mat1,3);
    W1 = mat1(:,:,1);
    %W1 = mat2gray(W1);

    img = imagesc(W1);
    colormap gray;
    clim(ptiles);

    pause(2);
    for K = 2:Q
        W1 = mat1(:,:,K);
        %W1 = mat2gray(W1);
        %W1 = imadjust(W1, [0,1],[]);
        set(img, 'CData', W1);
        %drawnow limitrate;
        drawnow();
        %caxis([0,1])
        pause(pause_time);
        %disp(K);
    end
end