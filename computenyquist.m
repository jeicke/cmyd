R = 60000;
phi = 30;
lamda = .03;
v = 50;
t = 2.5;
dy = [-v*t/2 v*t/2];
d = v/lamda .* cos(phi).*(dy./R)./sqrt(1+(dy./R).^2);%sin(atan2(dy,R))

min_prf = (max(d)-min(d)) * 2
