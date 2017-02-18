Z = load('assignmentShapeAnalysis/bone3D.mat');
Faces = Z.TriangleIndex;
Z = Z.shapesTotal;
% patch('Faces',Faces,'Vertices',Z(:,:,1)','FaceColor',rand(1,3));
[dim, n_points, n_pointsets] = size(Z);
T = mean(Z,2);
Z_trans = zeros(size(Z));
Z_normalised = zeros(size(Z));
for i=(1:n_pointsets)
    Z_trans(:,:,i) = Z(:,:,i) - T(:,:,i)*ones(1,n_points);
    Z_normalised(:,:,i) = Z_trans(:,:,i)/sqrt(sum(sum(Z_trans(:,:,i).^2)));
end
% y = Z_normalised(:,:,1);
% x = Z_normalised(:,:,2);
% z = aligned_pointset(y, x);
% 
% figure(1);
% % scatter(y(1,:),y(2,:),'b');
% % line(y(1,:),y(2,:),'Color','b');
% c = rand(1,3);
% patch('Faces',Faces,'Vertices',y','EdgeColor',c,'FaceColor',c,'FaceAlpha',.3);
% 
% hold on;
% c = rand(1,3);
% % scatter(x(1,:),x(2,:),'r');
% patch('Faces',Faces,'Vertices',x','EdgeColor',c,'FaceColor',c,'FaceAlpha',.3);
% hold on;
% c = rand(1,3);
% % scatter(z(1,:),z(2,:),'g');
% patch('Faces',Faces,'Vertices',z','EdgeColor',c,'FaceColor',c,'FaceAlpha',.3);

Z_aligned = Z_normalised;
z_mean = mean(Z_aligned, 3);
z_mean = z_mean - mean(z_mean,2)*ones(1,n_points);
z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));
for i= 1:10
    error = 0;
    for j = 1:n_pointsets
        Z_aligned(:,:,j) = aligned_pointset(z_mean, Z_aligned(:,:,j));
        error = error + sum(sum((z_mean-Z_aligned(:,:,j)).^2));
    end
    z_mean = mean(Z_aligned, 3);
    z_mean = z_mean - mean(z_mean,2)*ones(1,n_points);
    z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));
end

c = rand(1,3);
for i = 1:n_pointsets
    scatter3(Z(1,:,i),Z(2,:,i),Z(3,:,i),10,'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3));
    hold on;
end

patch('Faces',Faces,'Vertices',z_mean','EdgeColor','none','FaceColor','b','FaceAlpha',0.5);
hold off;
Z_aligned_reshape = zeros(dim*n_points,n_pointsets);
 for i = 1:n_pointsets
    Z_aligned_reshape(:,i) = reshape(Z_aligned(:,:,i),dim*n_points,1);%reshape(Z_aligned,dim*n_points,n_pointsets);
end
% z_mean = mean(Z_aligned, 3);
%z_mean_reshape = z_mean(:);

z_mean_reshape = mean(Z_aligned_reshape, 2);
Z_corr = Z_aligned_reshape - z_mean_reshape*ones(1,n_pointsets);
[V,D] = eig(Z_corr * Z_corr');

% z_mean = mean(Z_aligned, 3);
% z_mean = z_mean - mean(z_mean,2)*ones(1,n_points);
% z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));
% z_mean_reshape = z_mean(:);

figure();
A = z_mean_reshape;
A = reshape(A,dim,n_points);
patch('Faces',Faces,'Vertices',A','EdgeColor','none','FaceColor',[0,0,1],'FaceAlpha',0.4);
hold on;
B = z_mean_reshape + 0.5*sqrt(D(end,end))*V(:,end);
B = reshape(B,dim,n_points);
patch('Faces',Faces,'Vertices',B','EdgeColor','none','FaceColor',[0,0.1,0.9],'FaceAlpha',0.2);

hold on;
C = z_mean_reshape - 0.5*sqrt(D(end,end))*V(:,end);
C = reshape(C,dim,n_points);
patch('Faces',Faces,'Vertices',C','EdgeColor','none','FaceColor',[0.1,0,0.9],'FaceAlpha',0.2);


%modes = Z_corr*V;
