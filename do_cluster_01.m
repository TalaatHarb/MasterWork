function result = do_cluster_01(filename)
[v, f, ~, ~, ~] = stlread(filename);
[v, f]=patchslim(v, f);

% Number of vertices
num_points = size(v,1);

% Number of faces
num_faces = size(f,1);

% Array to store which points are in a cluster and which are not
PlaneID = zeros(num_points, 1);
currentPlane = 1;
tolerance = 9e-7;
distanceMeasure = zeros(num_points, num_points);
connectivity = zeros(num_points, num_points);
for faceCounter = 1:num_faces
    connectivity(f(faceCounter,1),f(faceCounter,2)) = 1;
    connectivity(f(faceCounter,2),f(faceCounter,1)) = 1;
    connectivity(f(faceCounter,2),f(faceCounter,3)) = 1;
    connectivity(f(faceCounter,3),f(faceCounter,2)) = 1;
    connectivity(f(faceCounter,3),f(faceCounter,1)) = 1;
    connectivity(f(faceCounter,1),f(faceCounter,3)) = 1;
end
for vertexCounter = 1:num_points
    % Ignoring already taken points
    if PlaneID(vertexCounter,1) ~= 0
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
        closestTwo = neighborhood(I(1:2,1),:);
        % Calculating the vertex plane
        normal = cross(closestTwo(1,:)-currentVertex, currentVertex-closestTwo(2,:));
        D = -dot(normal,currentVertex);
        for vertexCounter2 = 1:num_points
            % Ignoring already taken points
            if PlaneID(vertexCounter2,1) ~= 0
                continue;
            else
                % Add all verteces in the plane to the current cluster
                proximity_to_plane = abs((dot(normal, v(vertexCounter2,:))+ D)/norm(v(vertexCounter2,:)));
                if proximity_to_plane <= tolerance
                    PlaneID(vertexCounter2,1) = currentPlane;
                end
            end
        end
        currentPlane = currentPlane + 1;
    end
end
% Caculating the number of planes
p = max(PlaneID(:));
k = p;
ClusterID = PlaneID;
result = {num_points;...
    num_faces;...
    p;...
    k;...
    v;...
    f;...
    connectivity;...
    PlaneID;...
    ClusterID...
    };
for plane = 1:p
    result{9+plane} = v(PlaneID == plane,:);
end
end