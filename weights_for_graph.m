function [w, pmask, nmask] = weights_for_graph(Sw, prop,varargin)
%%% simple function to create linearly interpolated weights for graph
%%% plotting

vmax = varargin_parse(varargin,'vmax',6.25); %maximum feature size (width, or marker size)
vmin = varargin_parse(varargin,'vmin',0.25); %minimum feature size (width or marker size)

wmax = max(abs(Sw.(prop)));
wmin = min(abs(Sw.(prop)));
wrange = wmax - wmin;
w = (vmax-vmin)*(abs(Sw.(prop)) - wmin)./wrange + vmin;
pmask = Sw.(prop) > 0;
nmask = Sw.(prop) < 0;