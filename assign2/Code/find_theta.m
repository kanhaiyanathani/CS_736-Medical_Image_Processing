function [rrmse, theta] = find_theta(img)
for theta = (1:180)
    reconstructed_img= iradon(radon(img,theta+(0:149)),theta+(0:149));
    rrmse(theta) = RRMSE(img,reconstructed_img);
end
[~,theta] = min(rrmse);
end

