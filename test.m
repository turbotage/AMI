%%
rfmat = open('RF_MAT.mat').RF_MAT;

frame = 5;
Bmodes = sqrt(abs(hilbert(squeeze(rfmat(:,:,:)))));
shape = size(Bmodes)
Bmodes_f = reshape(Bmodes, shape(1)*shape(2), shape(3));

[U,S,V] = svd(Bmodes_f, 'econ');

Snew = S;
Snew(100:end, 100:end) = 0;
Snew(1:2,1:2) = 0;

Bmodes_fnew = U * Snew * V';
Bmodes_new = reshape(Bmodes_fnew, shape(1), shape(2), shape(3));

%figure(1);
%animate_stuff(Bmodes);
%figure(2);
%animate_stuff(Bmodes_new);
%%

disp('Hello World');
draw_pic(Bmodes, Bmodes_new);

function animate_stuff(mat)
    %movie(Bmodes, 1, 10);
    Q = size(mat,3);
    W = mat(:,:,1);
    h = pcolor(W);
    set(h, 'EdgeColor', 'none');
    drawnow();
    pause(0.3);
    for K = 2:Q
       W = mat(:,:,K);
       set(h, 'CData', W);
       drawnow limitrate;
       %pause(0.1);
       %disp(K);
    end
end

function draw_pic(mat1, mat2)
    figure(1);
    Q = size(mat1,3);
    W1 = mat1(:,:,1);
    W2 = mat2(:,:,2);
    subplot(1,2,1);
    img1 = imagesc(W1);
    colormap gray;
    subplot(1,2,2);
    img2 = imagesc(W2);
    colormap gray;
    
    for K = 2:Q
        W1 = mat1(:,:,K);
        W2 = mat2(:,:,K);
        set(img1, 'CData', W1);
        set(img2, 'CData', W2);
        %drawnow limitrate;
        drawnow();
        pause(0.02);
        disp(K);
    end
end
