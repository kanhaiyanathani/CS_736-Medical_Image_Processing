function [V,D] = modes(Z,z_mean,Z_aligned)
[dim,n_points,n_pointsets] = size(Z);
for i = 1:n_pointsets
    Z_aligned_reshape(:,i) = reshape(Z_aligned(:,:,i),dim*n_points,1);
end

z_mean_reshape = mean(Z_aligned_reshape, 2);
Z_corr = Z_aligned_reshape - z_mean_reshape*ones(1,n_pointsets);
[V,D] = eig(Z_corr * Z_corr');
D = D/300;
end