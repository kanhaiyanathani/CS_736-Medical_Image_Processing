function Rf = myRadonTrans(img,t,theta,s_step)
%t = -90:5:90;
%theta = 0:5:175;

Rf = zeros(length(t),length(theta));

for i = 1:length(t)
    for j = 1:length(theta)
        Rf(i,j) = myIntegration(img,t(i),theta(j),s_step);
    end
end

end
