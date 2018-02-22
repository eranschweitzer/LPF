function Gw = branchweights(Sb,varargin)

gmode     = varargin_parse(varargin,'GwAru','default');
bmode     = varargin_parse(varargin,'GwAmu','default');
bcmode    = varargin_parse(varargin,'bcmode','default');
bshmode   = varargin_parse(varargin,'bshmode','default');

Gw     = struct();
Gw.g   = genweights(abs(Sb.g),gmode);
Gw.b   = genweights(abs(Sb.b),bmode);
Gw.bc  = genweights(abs(Sb.bc),bcmode);
Gw.bsh = genweights(abs(Sb.bsh),bshmode);
