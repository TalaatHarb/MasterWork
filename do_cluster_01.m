function result = do_cluster_01(filename)
[v, f, ~, ~, ~] = stlread(filename);
[v, f]=patchslim(v, f);

% Number of vertices
num_points = size(v,1);

% Number of faces
num_faces = size(f,1);

% Array to store which points are in a cluster and which are not
clusterID = zeros(num_points, 1);
currentCluster = 1;
tolerance = 1e-6;
distanceMeasure = zeros(num_points, num_points);
for vertexCounter = 1:num_points
    % Ignoring already taken points
    if clusterID(vertexCounter,1) ~= 0
        continue;
    else
        currentVertex = v(vertexCounter, :);
        % Clustering
        tempList = [];
        for faceCounter = 1:num_faces
            % Getting which vertecies are connected to the current vertex
            if f(faceCounter, 1) == vertexCounter
                tempList = [tempList;f(faceCounter, 2);f(faceCounter, 3)];
            elseif f(faceCounter, 2) == vertexCounter
                tempList = [tempList;f(faceCounter, 1);f(faceCounter, 3)];
            elseif f(faceCounter, 3) == vertexCounter
                tempList = [tempList;f(faceCounter, 1);f(faceCounter, 2)];
            end
        end
        tempList = unique(tempList);
        num_connected = size(tempList, 1);
        neighborhood = v(tempList,:);
        distances = zeros(num_connected, 1);
        for neighbor = 1:num_connected
            % Calculating distance to neighbors
            if distanceMeasure(vertexCounter, tempList(neighbor)) == 0
                distances(neighbor,1) = norm (currentVertex - neighborhood(neighbor,:));
                distanceMeasure(vertexCounter, tempList(neighbor)) = distances(neighbor,1);
                distanceMeasure(tempList(neighbor), vertexCounter) = distances(neighbor,1);
            else
                distances(neighbor,1) = distanceMeasure(vertexCounter, tempList(neighbor));
            end
        end
        [~,I] = sort(distances(:,1));
        indexes_of_closest_2 = I(1:2,1);
        closestTwo = neighborhood(indexes_of_closest_2,:);
        % Calculating the vertex plane
        normal = cross(closestTwo(1,:)-currentVertex, currentVertex-closestTwo(2,:));
        D = -dot(normal,currentVertex);
        for vertexCounter2 = 1:num_points
            % Ignoring already taken points
            if clusterID(vertexCounter2,1) ~= 0
                continue;
            else
                % Add all verteces in the plane to the current cluster
                proximity_to_plane = abs((dot(normal, v(vertexCounter2,:))+ D)/norm(v(vertexCounter2,:)));
                if proximity_to_plane <= tolerance
                    clusterID(vertexCounter2,1) = currentCluster;
                end
            end
        end
        currentCluster = currentCluster + 1;
    end
end
% Caculating the number of cross sections
k = max(clusterID(:));

result = {num_points;num_faces;k;v;f;clusterID};
for cluster = 1:k
    result{6+cluster} = v(clusterID == cluster,:);
end
end