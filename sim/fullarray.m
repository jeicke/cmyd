function resp = fullarray(NESA,SBSA,spacing,SAsteer,FAsteer)
SAsteer = SAsteer*pi/180;
FAsteer = FAsteer*pi/180;
tot_el = 272;
% num el in subarray (NESA) SBSA (shift btwn sub arrays)
%NESA = 12;
%SBSA = 3;
% num of sub arrays (NSA)
NSA = 1 + floor((tot_el-NESA)/SBSA);
SABP = subarray(NESA,spacing,SAsteer);
FAwts = chebwin(NSA,40);
FAwts = FAwts/sum(FAwts); %normalize weights
SA_spacing = SBSA*spacing;
SA_posits = -((NSA - 1)/2)*SA_spacing:SA_spacing:((NSA-1)/2)*SA_spacing;
SA_posits = SA_posits';
theta = [-90:0.1:90]*pi/180;
SA_posits = SA_posits*2*pi;
e = exp(i*SA_posits*sin(theta));
FAwts = FAwts.*exp(i*SA_posits*sin(FAsteer));
resp = FAwts'*e;
resp = resp.*SABP;
resp = abs(resp).^2;
%figure
plot(theta,10*log10(resp));
hold on
axis([-1.75 1.75 -80 5]);
