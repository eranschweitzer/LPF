function phi = philimit(phi)

phimax = 0.5*(pi/2)^2;
phimin = 0;

phi(phi > phimax) = phimax;
phi(phi < phimin) = phimin;