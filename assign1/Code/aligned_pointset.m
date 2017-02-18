function [z_transformed] = aligned_pointset(ref,test) % test is pointset to be aligned with ref pointset
% ref and test (dim, n_points) are in preshape space i.e
% translated and normalized
dim = size(ref,1);

% Find Rotation matrix
cross_cov = test*ref';
[U,~,V] = svd(cross_cov);
R = V*U';
if det(R) == -1
    A = eye(dim);
    A(dim,dim) = -1;
    R = V*A*U';
end

z_transformed = R * test;  % test pointset aligned to the ref pointset - output
end