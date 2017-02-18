function [z_mean,Z_aligned] = mean_shape(Z)  % Z - given data set
[dim, n_points, n_pointsets] = size(Z);
T = mean(Z,2);
Z_trans = zeros(size(Z));
Z_normalised = zeros(size(Z));
for i=(1:n_pointsets)
    Z_trans(:,:,i) = Z(:,:,i) - T(:,:,i)*ones(1,n_points);  % Translation to origin
    Z_normalised(:,:,i) = Z_trans(:,:,i)/sqrt(sum(sum(Z_trans(:,:,i).^2))); % Normalisation
end

Z_aligned = Z_normalised;
z_mean = mean(Z_aligned, 3);   % initialisation of z_mean
z_mean = z_mean - mean(z_mean,2)*ones(1,n_points); % translation
z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));  % normalisation
for i= 1:10
    for j = 1:n_pointsets
        Z_aligned(:,:,j) = aligned_pointset(z_mean, Z_aligned(:,:,j));
    end
    z_mean = mean(Z_aligned, 3);
    z_mean = z_mean - mean(z_mean,2)*ones(1,n_points);
    z_mean = z_mean/sqrt(sum(sum(z_mean.^2)));
end
end