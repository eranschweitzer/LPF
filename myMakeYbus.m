function Y = myMakeYbus(F,T,Sb)

ys = Sb.g + 1i*Sb.b;

yff = (ys + 0.5*1i*Sb.bc)./Sb.tau.^2;
yft = -ys.*(exp(+1i*Sb.tshift)./Sb.tau);
ytf = -ys.*(exp(-1i*Sb.tshift)./Sb.tau);
ytt = ys + 0.5*1i*Sb.bc;

Y = F'*(diags(yff)*F + diags(yft)*T) + T'*(diags(ytt)*T + diags(ytf)*F) + diags(Sb.gsh + 1i*Sb.bsh);