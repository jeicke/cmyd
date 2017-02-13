v = 18;
f = 5E9;
c = 299792458;
lamda = c/f;
w = 1.2;
EL = 0;
theta = 0;
vangle = 90;
format  bank
df1 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta))
w = 2 * w;
df2 =  2 * v/lamda * cosd(EL) * 2 * (cosd(vangle-theta)-cosd(vangle + asind(lamda/w) - theta))

df2-df1

prf = 1500;
npulses = 256;

hz_spacing  = prf /npulses 