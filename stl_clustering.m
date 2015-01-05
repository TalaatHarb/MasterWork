%% Clearing Old Data
clear
clc
tic;
%% Opening the file
filename = 'multi.stl';
[v, f, n, c, stltitle] = stlread(filename);
[v, f]=patchslim(v, f);

% Number of vertices
num_points = size(v,1);

% Number of faces
num_faces = size(f,1);

% Display file parameters
disp(['File Name:                         ' filename]);
disp(['Number of vertices:                ' num2str(num_points)]);
disp(['Number of faces:                   ' num2str(num_faces)]);


%% Clustering
do_cluster_01;
%% Plotting the clusters

disp(['Number of cross sections:          ' num2str(k)]);

% http://www.mathworks.com/help/matlab/ref/colormap.html
colors = jet(k);

for cluster = 1:k
    X = v(clusterID == cluster,1);
    Y = v(clusterID == cluster,2);
    Z = v(clusterID == cluster,3);
    S = 16;
    C = repmat(colors(cluster,:),numel(X),1);
    scatter3(X,Y,Z,S,C);
    hold on;
end
patch('Faces',f,'Vertices',v);
alpha(0.5);
hold off;

toc;