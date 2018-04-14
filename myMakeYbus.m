function Y = myMakeYbus(F,T,Sb,varargin)

tronly = varargin_parse(varargin,'tronly',false);

ys = Sb.g + 1i*Sb.b;

if tronly
%     ytt = ys./Sb.tau;
%     yff = ytt;
    ytt = ys;
    yff = (ytt)./Sb.tau.^2;
else
    ytt = ys + 0.5*1i*Sb.bc;
    yff = (ytt)./Sb.tau.^2;
end
yft = -ys.*(exp(+1i*Sb.tshift)./Sb.tau);
ytf = -ys.*(exp(-1i*Sb.tshift)./Sb.tau);


Y = F'*(diags(yff)*F + diags(yft)*T) + T'*(diags(ytt)*T + diags(ytf)*F); 

if ~tronly
    Y = Y + diags(Sb.gsh + 1i*Sb.bsh);
end