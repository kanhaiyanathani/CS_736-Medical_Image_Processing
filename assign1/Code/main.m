%% Q1 Ellipses
Z = load('../Data/ellipses2D.mat');
Z = Z.pointSets;

[dim,n_points,n_pointsets] = size(Z);
% Q1 Part (a)
figure();
for i = 1:n_pointsets
    scatter(Z(1,:,i),Z(2,:,i),20,rand(1,3),'filled');
    hold on;
end
hold off;
title('Initial pointsets as given in dataset');

% Q1 Part (b)
[z_mean,Z_aligned] = mean_shape(Z); % Compute the mean shape
figure();
for i = 1:n_pointsets  % Plot aligned pointsets
    scatter(Z_aligned(1,:,i),Z_aligned(2,:,i),20,rand(1,3),'filled');
    hold on;
end
scatter(z_mean(1,:),z_mean(2,:),50,'r','filled');   % Plot computed mean shape
line(z_mean(1,:),z_mean(2,:),'Color','r','Linewidth',2);
title('Aligned pointsets and the mean shape');

% Q1 Part (c)
[V,D] = modes(Z,z_mean,Z_aligned); 
x = diag(D);
y = fliplr(x');  % Eigen values in decreasing order
figure();
plot(y);
title('Variances(eigenvalues;sorted)');

% Q1 Part (d)
z_mean = z_mean - mean(z_mean,2)*ones(1,n_points);
z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));
z_mean_reshape = z_mean(:);

figure();  % Plotting of mode 1
for i = 1:n_pointsets  % Plot aligned pointsets
    scatter(Z_aligned(1,:,i),Z_aligned(2,:,i),20,rand(1,3),'filled'); %de2bi(mod(i,8),3)
    hold on;
end
A = z_mean_reshape;    % Plot mean
A = reshape(A,dim,n_points);
scatter(A(1,:),A(2,:),'r','filled');
p1 = line(A(1,:),A(2,:),'Color','r','Linewidth',2.5);
hold on;
B = z_mean_reshape + 2*sqrt(D(end,end))*V(:,end);  
B = reshape(B,dim,n_points);
scatter(B(1,:),B(2,:),'g','filled');
p2 = line(B(1,:),B(2,:),'Color','g','Linewidth',2.5);
hold on;
C = z_mean_reshape - 2*sqrt(D(end,end))*V(:,end);
C = reshape(C,dim,n_points);
scatter(C(1,:),C(2,:),'b','filled');
p3 = line(C(1,:),C(2,:),'Color','b','Linewidth',2.5);
legend([p1 p2 p3],'Mean','Mean+2SD','Mean-2SD');
title('First principal mode of shape variation, mean and aligned pointset');


figure();  % Plotting of mode 2
for i = 1:n_pointsets  % Plot aligned pointsets
    scatter(Z_aligned(1,:,i),Z_aligned(2,:,i),20,rand(1,3),'filled'); %de2bi(mod(i,8),3)
    hold on;
end
A = z_mean_reshape;    % Plot mean
A = reshape(A,dim,n_points);
scatter(A(1,:),A(2,:),'r','filled');
p1 = line(A(1,:),A(2,:),'Color','r','Linewidth',2.5);
hold on;
B = z_mean_reshape + 2*sqrt(D(end-1,end-1))*V(:,end-1); % Mean + 2SD  
B = reshape(B,dim,n_points);
scatter(B(1,:),B(2,:),'g','filled');
p2 = line(B(1,:),B(2,:),'Color','g','Linewidth',2.5);
hold on;
C = z_mean_reshape - 2*sqrt(D(end-1,end-1))*V(:,end-1); % Mean -2SD
C = reshape(C,dim,n_points);
scatter(C(1,:),C(2,:),'b','filled');
p3 = line(C(1,:),C(2,:),'Color','b','Linewidth',2.5);
legend([p1 p2 p3],'Mean','Mean+2SD','Mean-2SD');
title('Second principal mode of shape variation, mean and aligned pointset');

%% Q2 Hands
Z = load('../Data/hands2D.mat');
Z = Z.shapes;

[dim,n_points,n_pointsets] = size(Z);
% Q2 Part (a)
figure();
for i = 1:n_pointsets
    scatter(Z(1,:,i),Z(2,:,i),20,rand(1,3),'filled');
    hold on;
end
hold off;
title('Initial pointsets as given in dataset');

% Q2 Part (b)
[z_mean,Z_aligned] = mean_shape(Z); % Compute the mean shape
figure();
for i = 1:n_pointsets  % Plot aligned pointsets
    scatter(Z_aligned(1,:,i),Z_aligned(2,:,i),20,rand(1,3),'filled'); %de2bi(mod(i,8),3)
    hold on;
end
scatter(z_mean(1,:),z_mean(2,:),50,'r','filled');   % Plot computed mean shape
line(z_mean(1,:),z_mean(2,:),'Color','r','Linewidth',2);
title('Aligned pointsets and the mean shape');

% Q2 Part (c)
[V,D] = modes(Z,z_mean,Z_aligned);
x = diag(D);
y = fliplr(x');
figure();
plot(y);
title('Variances(eigenvalues;sorted)');

% Q2 Part (d)
z_mean = z_mean - mean(z_mean,2)*ones(1,n_points);
z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));
z_mean_reshape = z_mean(:);

figure();  % Plotting of mode 1
for i = 1:n_pointsets  % Plot aligned pointsets
    scatter(Z_aligned(1,:,i),Z_aligned(2,:,i),20,rand(1,3),'filled'); %de2bi(mod(i,8),3)
    hold on;
end
A = z_mean_reshape;    % Plot mean
A = reshape(A,dim,n_points);
scatter(A(1,:),A(2,:),'r','filled');
p1 = line(A(1,:),A(2,:),'Color','r','Linewidth',2.5);
hold on;
B = z_mean_reshape + 2*sqrt(D(end,end))*V(:,end);  % Mean + 2SD
B = reshape(B,dim,n_points);
scatter(B(1,:),B(2,:),'g','filled');
p2 = line(B(1,:),B(2,:),'Color','g','Linewidth',2.5);
hold on;
C = z_mean_reshape - 2*sqrt(D(end,end))*V(:,end);
C = reshape(C,dim,n_points);
scatter(C(1,:),C(2,:),'b','filled');
p3 = line(C(1,:),C(2,:),'Color','b','Linewidth',2.5);  % Mean - 2SD
legend([p1 p2 p3],'Mean','Mean+2SD','Mean-2SD');
title('First principal modes of shape variation, mean and aligned pointset');

figure();  % Plotting of mode 1
for i = 1:n_pointsets  % Plot aligned pointsets
    scatter(Z_aligned(1,:,i),Z_aligned(2,:,i),20,rand(1,3),'filled'); %de2bi(mod(i,8),3)
    hold on;
end
A = z_mean_reshape;    % Plot mean
A = reshape(A,dim,n_points);
scatter(A(1,:),A(2,:),'r','filled');
p1 = line(A(1,:),A(2,:),'Color','r','Linewidth',2.5);
hold on;
B = z_mean_reshape + 2*sqrt(D(end-1,end-1))*V(:,end-1); % Mean + 2SD  
B = reshape(B,dim,n_points);
scatter(B(1,:),B(2,:),'g','filled');
p2 = line(B(1,:),B(2,:),'Color','g','Linewidth',2.5);
hold on;
C = z_mean_reshape - 2*sqrt(D(end-1,end-1))*V(:,end-1);  % Mean - 2SD
C = reshape(C,dim,n_points);
scatter(C(1,:),C(2,:),'b','filled');
p3 = line(C(1,:),C(2,:),'Color','b','Linewidth',2.5);
legend([p1 p2 p3],'Mean','Mean+2SD','Mean-2SD');
title('Second principal mode of shape variation, mean and aligned pointset');


%% Q3 Bone
Z = load('../Data/bone3D.mat');
Faces = Z.TriangleIndex;
Z = Z.shapesTotal;
[dim,n_points,n_pointsets] = size(Z);

% Q3 Part (a)
set(figure(), 'Position', [1 1 800 500 ]);
for i = 1:n_pointsets
    scatter3(Z(1,:,i),Z(2,:,i),Z(3,:,i),15,'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3));
    hold on;
end
hold off;
title('Initial pointsets as given in dataset');

% Q3 Part (b)
[z_mean,Z_aligned] = mean_shape(Z); % Compute the mean shape
set(figure(), 'Position', [1 1 800 500 ]);
for i = 1:n_pointsets  % Plot aligned pointsets
    scatter3(Z_aligned(1,:,i),Z_aligned(2,:,i),Z_aligned(3,:,i),15,'MarkerEdgeColor','none','MarkerFaceColor',rand(1,3));
    hold on;
end
% Plot computed mean shape
patch('Faces',Faces,'Vertices',z_mean','EdgeColor','none','FaceColor','b','FaceAlpha',0.5);
title('Aligned pointsets and the mean shape');

% Q3 Part (c)
[V,D] = modes(Z,z_mean,Z_aligned);
x = diag(D);
y = fliplr(x');
figure();
plot(y);
title('Variances(eigenvalues;sorted)');

% Q3 Part (d)
z_mean = z_mean - mean(z_mean,2)*ones(1,n_points);
z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));
z_mean_reshape = z_mean(:);

set(figure(), 'Position', [1 1 1000 500 ]);
B = z_mean_reshape + 2*sqrt(D(end,end))*V(:,end);  
B = reshape(B,dim,n_points);
p2 = patch('Faces',Faces,'Vertices',B','EdgeColor','k','FaceColor','b','FaceAlpha',0.3);
title('Mean+2SD');

set(figure(), 'Position', [1 1 1000 500 ]);
A = z_mean_reshape;    % Plot mean
A = reshape(A,dim,n_points);
p1 = patch('Faces',Faces,'Vertices',A','EdgeColor','k','FaceColor','b','FaceAlpha',0.4);
title('Mean shape');

set(figure(), 'Position', [1 1 1000 500 ]);
C = z_mean_reshape - 2*sqrt(D(end,end))*V(:,end);
C = reshape(C,dim,n_points);
p3 = patch('Faces',Faces,'Vertices',C','EdgeColor','k','FaceColor','b','FaceAlpha',0.3);
title('Mean-2SD');

