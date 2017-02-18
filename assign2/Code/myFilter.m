function [filtered_h] = myFilter(h, filter, L, n_fft)

Fh = fft(h,n_fft);
% w = [0:n_fft/2-1, n_fft/2-1:0];
w = -n_fft/2:n_fft/2-1;
rect = rectangularPulse(-L,L-1,w)>0;

if strcmp(filter,'Ram-Lak')
    Cw = 1;
elseif strcmp(filter,'Shepp-Logan')
    Cw = sinc(0.5*w/L);
else
    Cw = cos(0.5*pi*w/L);
end

Aw = abs(w).*Cw.*rect*2/n_fft; 
Aw = [Aw(n_fft/2+1:end), Aw(1:n_fft/2)];
% plot(Aw);
Aw_m = Aw'*ones(1,size(h,2));

filtered_h = real(ifft(Aw_m.*Fh));
filtered_h(length(h)+1:end,:)=[];
end
