%% Clearing Old Data
clear all;
clf;
clc;
tic;
%% Opening the file
filename = 'test.stl';

%% Clustering
result = do_cluster_02(filename);
return
num_points = result{1};
num_faces = result{2};
p = result{3};
k = result{3};
v = result{5};
f = result{6};
connectivity = result{7};
PlaneID = result{8};
ClusterID = result{9};

%% Plotting the clusters
% Display file parameters
disp(['File Name:                         ' filename]);
disp(['Number of vertices:                ' num2str(num_points)]);
disp(['Number of faces:                   ' num2str(num_faces)]);
disp(['Number of cross sections:          ' num2str(p)]);

% http://www.mathworks.com/help/matlab/ref/colormap.html
colors = jet(k);

for cluster = 1:k
    
    current_cluster = result{9+cluster};
    
    % Displaying some information about the cluster
    fprintf('\n');
    disp(['Cluster ' num2str(cluster)]);
    disp(['Cluster Center:                    ' num2str(mean(current_cluster))]);
    disp(['Cluster Radius:                    ' num2str(norm(std(current_cluster)))]);
    X = current_cluster(:,1);
    Y = current_cluster(:,2);
    Z = current_cluster(:,3);
    S = 16;
    C = repmat(colors(cluster,:),numel(X),1);
    
    scatter3(X,Y,Z,S,C);
    hold on;
    
    clearvars X Y Z S C
end

clearvars colors cluster current_cluster

patch('Faces',f,'Vertices',v);
alpha(0.5);
hold off;

fprintf('\n');

toc;