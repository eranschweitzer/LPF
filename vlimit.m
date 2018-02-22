function [u0,busmask] = vlimit(u0,E,varargin)

vmax = varargin_parse(varargin,'vmax',1.3);
vmin = varargin_parse(varargin,'vmin',0.5);
uopt = varargin_parse(varargin,'uopt','flat');
% u0(u0 > log(vmax)) = 0;
% u0(u0 < log(vmin)) = 0;

if ~strcmp(uopt,'none')
    busmask = (u0 > log(vmax)) | (u0 < log(vmin));
    if any(busmask)
        switch uopt
            case 'flat'
                uavg = mean(u0(~busmask))*ones(sum(busmask),1);
            case 'avg'
                uavg = ufix(u0,E,busmask);
            case 'null'
                uavg = zeros(sum(busmask),1);
        end
        u0(busmask) = uavg;
    end
else
    busmask = false;
end