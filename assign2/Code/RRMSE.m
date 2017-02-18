function [out] = RRMSE(img1,img2)
out = norm(img1-img2(2:end-1,2:end-1))/norm(img1);
end