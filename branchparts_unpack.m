function [g,b,bc,tau,tshift,gsh,bsh] = branchparts_unpack(S,btype)
if nargin == 1
    btype = 0;
end
g      = S.g;
bc     = S.bc;
tau    = S.tau;
tshift = S.tshift;
gsh    = S.gsh;
bsh    = S.bsh;

if btype == 0
    b  = S.b;
elseif btype == 1
    b  = S.balt;
else
    error('Inocorrect btype')
end