function draw_pic2(mat1, mat2, delay, climc)
% plots 2 images with specified delay and option of including
% figure number

if nargin < 4
    climc = [0 1];
end

lowi = 0.0;
highi = 0.1;

figure
if size(mat1,3)==3
    Q = size(mat1,4);
    W1 = mat1(:,:,:,1);
else
    Q = size(mat1,3);
    W1 = mat1(:,:,1);
    W1 = mat2gray(W1);
end
W2 = mat2(:,:,:,1);

% Maximize window.
g = gcf;
% 		g.WindowState = 'maximized';
drawnow;

subplot(1,2,1);
if size(mat1,3)==3
    img1 = imagesc(W1, [lowi, highi]);
    title('At rest')
else
    img1 = imagesc(W1, [lowi, highi]);colormap gray;
    title('B-mode')
end

subplot(1,2,2);
img2 = imagesc(W2, [lowi, highi]);
if size(mat1,3)==3
    title('Contracted')
    sgtitle(['Comparison of activity between fully contracted and ' ...
        'at rest'])
else
    title('H-scan')
end

pause(2);

for K = 2:Q
    if size(mat1,3)==3
        W1 = mat1(:,:,:,K);
    else
        W1 = mat1(:,:,K);
        W1 = mat2gray(W1);
    end
    W1 = imadjust(W1, [0,1],[]);

    W2 = mat2(:,:,:,K);
    W2 = imadjust(W2, [0,1],[]);

    set(img1, 'CData', W1);
    set(img2, 'CData', W2);

    %drawnow limitrate;
    drawnow();
    clim(climc)
    pause(delay);
    disp(K);
end
end