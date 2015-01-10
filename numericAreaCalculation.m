function [center, radius, area, numericArea] = numericAreaCalculation(cluster)
% Calculates the area inside a convex curve assuming the area can be
% divided into trinagles from the center of the cluster of point composing
% the curve

center = mean(cluster);
radius = norm(std(cluster));
area = pi * radius * radius;



clusterSize = size(cluster,1);

numericArea = 0;
for i = 1:(clusterSize - 1)
    vector1 = cluster(i, :) - center;
    vector2 = cluster(i+1, :) - center;
    
    triArea = 0.5 * norm(cross(vector1,vector2));
    
    numericArea = numericArea + triArea;
    
end
vector1 = cluster(end, :) - center;
vector2 = cluster(1, :) - center;

triArea = 0.5 * norm(cross(vector1,vector2));

numericArea = numericArea + triArea;

end