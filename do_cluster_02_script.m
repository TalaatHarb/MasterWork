clearvars -except filename
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
dissimilarity = sparse(num_points, num_points);
connectivity = sparse(num_points, num_points);
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
        [~,tempList] = find(connectivity(vertexCounter,:)>0);
        tempList = tempList.';
        num_connected = size(tempList, 1);
        neighborhood = v(tempList,:);
        distances = zeros(num_connected, 1);
        for neighbor = 1:num_connected
            % Calculating distance to neighbors
            if dissimilarity(vertexCounter, tempList(neighbor)) == 0
                distances(neighbor,1) = norm (currentVertex -...
                    neighborhood(neighbor,:));
                dissimilarity(vertexCounter, tempList(neighbor)) =...
                    distances(neighbor,1);
                dissimilarity(tempList(neighbor), vertexCounter) =...
                    distances(neighbor,1);
            else
                distances(neighbor,1) = dissimilarity(vertexCounter,...
                    tempList(neighbor));
            end
        end
        [~,I] = sort(distances(:,1));
        closestTwo = neighborhood(I(1:2,1),:);
        % Calculating the vertex plane
        normal = cross(closestTwo(1,:)-currentVertex, currentVertex -...
            closestTwo(2,:));
        D = -dot(normal,currentVertex);
        for vertexCounter2 = 1:num_points
            % Ignoring already taken points
            if PlaneID(vertexCounter2,1) ~= 0
                continue;
            else
                % Add all verteces in the plane to the current cluster
                proximity = abs((dot(normal, v(vertexCounter2,:))+...
                    D)/norm(v(vertexCounter2,:)));
                if proximity <= tolerance
                    PlaneID(vertexCounter2,1) = currentPlane;
                end
            end
        end
        currentPlane = currentPlane + 1;
    end
end

% Caculating the number of planes
p = max(PlaneID(:));
planes = cell(p,1);

ClusterID = zeros(num_points, 1);
currentCluster = 1;
for plane = 1:p
    
    planes{plane,1} = v(PlaneID == plane,:);
    currentPlane = planes{plane,1};
    currentVertices = find(PlaneID(:,1) == plane);
    
    nextVertex = 1;
    while(~all(ClusterID(currentVertices,1) ~= 0))
        
        % Ignoring already taken points
        currentVertexID = currentVertices(nextVertex,1);
        if (ClusterID(currentVertexID,1) ~= 0) || isempty(nextVertex)
            currentCluster = currentCluster + 1;
            [~,nextVertex] = min(ClusterID(currentVertices,:));
            continue;
        else
            currentVertex = currentPlane(nextVertex, :);
            % Clustering
            ClusterID(currentVertexID,1) = currentCluster;
            
            [~,tempList] = find(connectivity(currentVertexID,:)>0);
            tempList = tempList(ClusterID(tempList) == 0).';
            
            if isempty(tempList)
                currentVertexID = 1;
                continue;
            end
            
            num_connected = size(tempList, 1);
            neighborhood = v(tempList,:);
            distances = zeros(num_connected, 1);
            for neighbor = 1:num_connected
                % Calculating distance to neighbors
                if dissimilarity(currentVertexID, tempList(neighbor)) == 0
                    distances(neighbor,1) = norm (currentVertex -...
                        neighborhood(neighbor,:));
                    dissimilarity(currentVertexID, tempList(neighbor)) =...
                        distances(neighbor,1);
                    dissimilarity(tempList(neighbor), currentVertexID) =...
                        distances(neighbor,1);
                else
                    distances(neighbor,1) = dissimilarity(currentVertexID...
                        , tempList(neighbor));
                end
            end
            [~,I] = min(distances(:,1));
            closestVertex = tempList(I,1);
        end
        nextVertex = find(currentVertices == closestVertex);
    end
    currentCluster = currentCluster + 1;
end

k = max(ClusterID(:));
clusters = cell(k,1);
for cluster = 1:k
    clusters{cluster,1} = v(ClusterID == cluster,:);
end

result = {num_points;...
    num_faces;...
    p;...
    k;...
    v;...
    f;...
    connectivity;...
    PlaneID;...
    ClusterID;...
    planes;...
    clusters;...
    dissimilarity...
    };
clearvars -except result filename