function Outliner(matrix,threshold_l,threshold_u,color_box)

[nRows, nCols] = size(matrix);

for r = 1:nRows
    for c = 1:nCols
        
        if matrix(r,c) > threshold_l && matrix(r,c) < threshold_u
            
            % Top edge if no neighbor above or out of bounds
            if r == 1 || matrix(r-1,c) <= threshold_l || matrix(r-1,c) >= threshold_u
                line([c-0.5, c+0.5], [r-0.5, r-0.5], ...
                     'Color',color_box,'LineWidth',2);
            end
            
            % Bottom edge if no neighbor below or out of bounds
            if r == nRows || matrix(r+1,c) <= threshold_l || matrix(r+1,c) >= threshold_u
                line([c-0.5, c+0.5], [r+0.5, r+0.5], ...
                     'Color',color_box,'LineWidth',2);
            end
            
            % Left edge if no neighbor to the left or out of bounds
            if c == 1 || matrix(r,c-1) <= threshold_l || matrix(r,c-1) >= threshold_u
                line([c-0.5, c-0.5], [r-0.5, r+0.5], ...
                     'Color',color_box,'LineWidth',2);
            end
            
            % Right edge if no neighbor to the right or out of bounds
            if c == nCols || matrix(r,c+1) <= threshold_l || matrix(r,c+1) >= threshold_u
                line([c+0.5, c+0.5], [r-0.5, r+0.5], ...
                     'Color',color_box,'LineWidth',2);
            end
            
        end
    end
end

end