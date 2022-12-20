function draw_pic2(mat1, mat2, delay, fignum)
% plots 2 images with specified delay and option of including
% figure number
	lowi = 0.0;
	highi = 0.1;

	if fignum
		figure(fignum);
	else 
		figure
	end
		Q = size(mat1,3);
		
	if size(mat1,3)==3
		W1 = mat1(:,:,:,1);
	else
		W1 = mat1(:,:,1);
		W1 = mat2gray(W1);
	end
    	W2 = mat2(:,:,:,1);

		% Maximize window.
		g = gcf;
% 		g.WindowState = 'maximized';
		drawnow;
	
    	subplot(1,2,1);
    	img1 = imagesc(W1, [lowi, highi]);colormap gray;
        if size(mat1,3)==3
		    title('At rest')
        else
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
        clim([6.5e-3,0.65])
        pause(delay);
        disp(K);
		end
end