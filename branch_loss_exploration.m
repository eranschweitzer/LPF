% branch losses

vmax = 1.1; % p.u.
vmin = 0.9; % p.u.
delta_max = 45; % deg

v1 = vmin:0.01:vmax;
y = 1/(1 + 1i*10);
for v = v1
    v2mag = v1;
    theta = (-delta_max:delta_max)*pi/180;
    [V,T] = meshgrid(v2mag,theta);
    Sloss = y*(v - V.*exp(1i*T)).*conj(v - V.*exp(1i*T));   
end