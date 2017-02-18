function [integral] = myIntegration(img, t, theta,s_step)
[r,c] = size(img);
[X,Y] = meshgrid(-c/2+1:c/2, r/2:-1:-r/2+1);
s = -(sqrt(r^2+c^2))/2+1:s_step:(sqrt(r^2+c^2))/2;
Xq = t*cos(theta*pi/180) - s* sin(theta*pi/180);
Yq = t*sin(theta*pi/180) + s* cos(theta*pi/180);

output = interp2(X,Y,img,Xq,Yq);
output(isnan(output)) = 0;

integral = sum(output);
end

