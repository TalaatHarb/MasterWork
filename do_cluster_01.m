% Array to store which points are in a cluster and which are not
cluster = zeros(num_points, 1);
currentCluster = 1;
tolerance = 1e-6;;
for vertexCounter = 1:num_points
    % Ignoring already taken points
    if cluster(vertexCounter,1) ~= 0
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
            distances(neighbor,1) = norm (currentVertex - neighborhood(neighbor,:));
        end
        [M,I] = sort(distances(:,1));
        indexes_of_closest_2 = I(1:2,1);
        closestTwo = neighborhood(indexes_of_closest_2,:);
        % Calculating the vertex plane
        normal = cross(closestTwo(1,:)-currentVertex, currentVertex-closestTwo(2,:));
        D = -dot(normal,currentVertex);
        for vertexCounter2 = 1:num_points
            % Ignoring already taken points
            if cluster(vertexCounter2,1) ~= 0
                continue;
            else
                % Add all verteces in the plane to the current cluster
                proximity_to_plane = dot(normal, v(vertexCounter2,:))+ D;
                if proximity_to_plane <= tolerance
                    cluster(vertexCounter2,1) = currentCluster;
                end
            end
        end
        currentCluster = currentCluster + 1;
    end
end
% Caculating the number of clusters
k = max(cluster(:))