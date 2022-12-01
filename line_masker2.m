function masked = line_masker2(mask1, mask2)
    xmax = size(mask1,1);
    ymax = size(mask1,2);
    masked = zeros(xmax, ymax);

    for i=1:ymax
        if mask1(1,i) < 0.5
            continue;
        end
        
        for j=1:xmax
           if mask2(j,i) > 0.5
                masked(j,i) = 1;
           else
               break;
           end
        end

        if mask1(xmax,i) < 0.5
            continue;
        end

        for j=flip(1:xmax)
           if mask2(j,i) > 0.5
                masked(j,i) = 1;
           else
               break;
           end
        end
    end
end