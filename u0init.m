function u0 = u0init(vg,bidx,varargin)

udefault = varargin_parse(varargin,'udefault','uref');
if iscell(udefault)
    u0 = cell(length(udefault),1);
    for k = 1:length(udefault)
        u0{k} = u0init(vg,bidx,'udefault',udefault{k});
    end
    return
end

% initialize the u0 vector using the voltage set point vector vg
N = length(vg);
upv = log(vg(bidx.pv));
uref = log(vg(bidx.ref));

if strcmp(udefault,'uref')
    udefault = uref;
end

u0 = udefault*ones(N,1);
u0(bidx.pv)  = upv;
u0(bidx.ref) = uref;