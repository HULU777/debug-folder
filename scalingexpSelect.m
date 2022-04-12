function [scaled,sortPselected,parents] = scalingexpSelect(fitness,c)
    parentSize = size(fitness,1);
    avgf = mean(fitness);
    maxf = max(fitness);
    scaled = exp(c*(fitness-avgf)/(maxf-avgf));
    [parents,sortPselected] = SimultaneousscaleSelect(scaled,parentSize);
end