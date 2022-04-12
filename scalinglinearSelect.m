function [scaled,sortPselected,parents] = scalinglinearSelect(fitness,c)
    minf = min(fitness);
    avgf = mean(fitness);
    maxf = max(fitness);
    threshold = (c*avgf-maxf)/(c-1);
    if minf > threshold
        alpha = (c-1)*avgf/(maxf-avgf);
        beta = (maxf-c*avgf)*maxf/(maxf-avgf);
    else
        alpha = avgf/(maxf-avgf);
        beta = -minf*avgf/(avgf-minf);
    end
    scaled = alpha*fitness+beta; 
    
    parentSize = size(fitness,1);
    [parents,sortPselected] = SimultaneousscaleSelect(scaled,parentSize);
end
