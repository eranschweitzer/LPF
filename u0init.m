function u0 = u0init(vg,bidx,varargin)

udefault = varargin_parse(varargin,'udefault','uref');
unearest = varargin_parse(varargin,'unearest',[]);
if strcmp(udefault,'nearest') && isempty(unearest)
    error('u0init: when using the nearest pv/ref bus method a vector of the nearest bus unearest must be provided.')
end
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

if strcmp(udefault,'nearest')
    u0 = ones(N,1);
else
    u0 = udefault*ones(N,1);
end
u0(bidx.pv)  = upv;
u0(bidx.ref) = uref;

if strcmp(udefault,'nearest')
    upq = u0(bidx.pq);
    upq(unearest == 0) = uref;
    upq(unearest ~= 0) = upv(unearest(unearest ~= 0));
    u0(bidx.pq) = upq;
end