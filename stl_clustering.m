%% Clearing Old Data
clear
clc

%% Opening the file
filename = 'cyl.stl';
[v, f, n, c, stltitle] = stlread(filename);
[v, f]=patchslim(v, f);

% Number of vertices
num_points = size(v,1);

% Number of faces
num_faces = size(f,1);

% Display file parameters
disp(['File Name:          ' filename]);
disp(['Number of vertices: ' num2str(num_points)]);
disp(['Number of faces:    ' num2str(num_faces)]);


%% Clustering
do_cluster_01;
return
%% Plotting the clusters
cmap = hsv(k);
for cluster = 1:k
    plot3(m(cluster,1),m(cluster,2),m(cluster,3),'o','Color',cmap(cluster,:));
    hold on;
    plot3(v(id(1:count(cluster),cluster),1),v(id(1:count(cluster),cluster),2),v(id(1:count(cluster),cluster),3),'+','Color',cmap(cluster,:));
end
patch('Faces',f,'Vertices',v,'FaceVertexCData',c);
hold off;