%% Clearing Old Data
clear all
close all
clc
tic;
%% Opening the file
filename = 'multi.stl';

%% Clustering
result = do_cluster_01(filename);

num_points = result{1};
num_faces = result{2};
k = result{3};
v = result{4};
f = result{5};
clusterID = result{6};

%% Plotting the clusters
% Display file parameters
disp(['File Name:                         ' filename]);
disp(['Number of vertices:                ' num2str(num_points)]);
disp(['Number of faces:                   ' num2str(num_faces)]);
disp(['Number of cross sections:          ' num2str(k)]);

% http://www.mathworks.com/help/matlab/ref/colormap.html
colors = jet(k);

for cluster = 1:k
    current_cluster = result{6+cluster};
    X = current_cluster(:,1);
    Y = current_cluster(:,2);
    Z = current_cluster(:,3);
    S = 16;
    C = repmat(colors(cluster,:),numel(X),1);
    scatter3(X,Y,Z,S,C);
    hold on;
    clearvars X Y Z S C
end
clearvars colors
patch('Faces',f,'Vertices',v);
alpha(0.5);
hold off;

toc;