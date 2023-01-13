function draw_pic(mat1, mat2, mat3, mat4, delay, climc)
% Rebeccas adaption of draw_pic for H-scan. Plots 3 gray or 4 rgb images with specified delay and option of including
% color limits (has default.)

if nargin < 6
    climc = [0 1];
end

lowi = 0.0;
highi = 0.1;

figure
% rgb or not condition
if size(mat2,3)==3
    Q = size(mat2,4);

    W1 = mat1(:,:,:,1);

    W2 = mat2(:,:,:,1);

    W3 = mat3(:,:,:,1);

    W4 = mat4(:,:,:,1);

    % Maximize window.
    g = gcf;
    % g.WindowState = 'maximized';
    drawnow;

    subplot(2,2,1);
    img1 = imagesc(W1, [lowi, highi]);
    title('Original Bmode color channel')

    subplot(2,2,2);
    img2 = imagesc(W2, [lowi, highi]);
    title('Low pass filtered Bmodes color channel')

    subplot(2,2,3);
    img3 = imagesc(W3, [lowi, highi]);
    title('High pass filtered Bmodes color channel')

    subplot(2,2,4);
    img4 = imagesc(W4, [lowi, highi]);
    title('Rgb coded bandpass filtered Bmode')
    sgtitle('H-scan images compared to the original B-modes');

    pause(2);

    for K = 2:Q
        W1 = mat1(:,:,:,K);
        W1 = imadjust(W1, [0,1],[]);

        W2 = mat2(:,:,:,K);
        W2 = imadjust(W2, [0,1],[]);

        W3 = mat3(:,:,:,K);
        W3 = imadjust(W3, [0,1],[]);

        W4 = mat4(:,:,:,K);
        W4 = imadjust(W4, [0,1],[]);

        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
        set(img3, 'CData', W3);
        set(img4, 'CData', W4);

        %drawnow limitrate;
        drawnow();
        clim(climc)
        pause(delay);
        disp(K);
    end

else
    Q = size(mat1,3);

    W1 = mat1(:,:,1);
    W1 = mat2gray(W1);

    W2 = mat2(:,:,1);
    W2 = mat2gray(W2);

    W3 = mat3(:,:,1);
    W3 = mat2gray(W3);

    % Maximize window.
    g = gcf;
    %     	g.WindowState = 'maximized';
    drawnow;

    subplot(1,3,1);
    img1 = imagesc(W1, [lowi, highi]);colormap gray;
    title('Original Bmode')

    subplot(1,3,2);
    img2 = imagesc(W2, [lowi, highi]);
    title('Low pass filtered Bmode')

    subplot(1,3,3);
    img3 = imagesc(W3, [lowi, highi]);
    title('High pass filtered Bmode')

    sgtitle('H-scan images compared to the original B-modes');

    pause(2);

    for K = 2:Q
        W1 = mat1(:,:,K);
        W1 = mat2gray(W1);
        W1 = imadjust(W1, [0,1],[]);

        W2 = mat2(:,:,K);
        W2 = mat2gray(W2);
        W2 = imadjust(W2, [0,1],[]);

        W3 = mat3(:,:,K);
        W3 = mat2gray(W3);
        W3 = imadjust(W3, [0,1],[]);

        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
        set(img3, 'CData', W3);

        %drawnow limitrate;
        drawnow();
        clim([0,1])
        pause(delay);
        disp(K);
    end
end
end