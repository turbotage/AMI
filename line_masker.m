function masked = line_masker(mask, begin)
    fillLength = 20;
    fillRatio = 0.9;
    perc = 0.05;
    masked = mask;
    xmax = size(masked,1);
    u = xmax*(1-perc);
    l = xmax*perc;

    [xi,yi] = find(masked);

    len = length(xi);
    halflen = round(len/2);
    for i=1:halflen
        xpos = xi(i);
        ypos = yi(i);
        

        if begin && ((xpos<u)&&(xpos>l))
            masked(xpos,ypos) = 0;
            continue;
        end

        clampedRange = max(1,xpos-fillLength):min(xmax,xpos+fillLength);
        clampedRangeValues = masked(clampedRange,ypos);
        
        count = sum(clampedRangeValues);
        indices = find(clampedRangeValues);
        firstIndex = min(indices);
        lastIndex = max(indices);
        
        %fprintf("count=%d, firstIndex=%d, lastIndex=%d", count, firstIndex, lastIndex);
    
        clampedRange2 = max(1,xpos-fillLength+firstIndex):min(xmax,xpos+lastIndex);
        if count / fillLength > fillRatio
            masked(clampedRange2,ypos) = 1;
        else
            %masked(clampedRange,ypos) = 0;
        end
    end
    
    for i=flip(max(halflen,1):len)
        xpos = 0;
        ypos = 0;
        try
            xpos = xi(i);
            ypos = yi(i);
        catch
            fprintf('i=%d', i);
        end

        if begin && ((xpos<u)&&(xpos>l))
            masked(xpos,ypos) = 0;
            continue;
        end

        clampedRange = max(1,xpos-fillLength):min(xmax,xpos+fillLength);
        clampedRangeValues = masked(clampedRange,ypos);
        
        count = sum(clampedRangeValues);
        indices = find(clampedRangeValues);
        firstIndex = min(indices);
        lastIndex = max(indices);
        
        %fprintf("count=%d, firstIndex=%d, lastIndex=%d", count, firstIndex, lastIndex);
    
        clampedRange2 = max(1,xpos-fillLength+firstIndex):min(xmax,xpos+lastIndex);
        if count / fillLength > fillRatio
            masked(clampedRange2,ypos) = 1;
        else
            %masked(clampedRange,ypos) = 0;
        end
    end

end