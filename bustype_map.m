function bidx = bustype_map(mpc,bus_with_ongen,gstat,nmap,varargin)
%returns boolean maps indicating which bus is of which type
%bus with on gen is an Nx1 boolean vector with 1 if the bus has a generator that
%is on and 0 otherwise.

newref = varargin_parse(varargin,'newref',1);

define_constants;
bidx = struct();
bidx.pq  = mpc.bus(:,BUS_TYPE) == PQ;
bidx.pv  = mpc.bus(:,BUS_TYPE) == PV;
bidx.ref = mpc.bus(:,BUS_TYPE) == REF;

% if a bus is specified as PV but has no generator attached change it to PQ
% same for ref bus (though this really shouldn't happen normally)
pv2pq  = find(bidx.pv  & ~bus_with_ongen);
bidx.pv(pv2pq)   = false;
bidx.pq(pv2pq)   = true;

ref2pq = find(bidx.ref & ~bus_with_ongen);
bidx.ref(ref2pq) = false;
bidx.pq(ref2pq) = true;

if ~any(bidx.ref)
    warning(['No Reference Bus: ',...
        'Possibly designated reference bus had no generator.\n',...
        'Reassigning using option %d'], newref)
    
    switch newref
        case 1
            %option 1: like matpower select first pv bus
            ref = find(bidx.pv,1);
        case 2
            %option 2: select dispatched generator with most headroom
            [~,idx] = max((mpc.gen(:,PMAX) - mpc.gen(:,PG)).*gstat);
            ref = full(nmap(mpc.gen(idx,1)));
        case 3 
            %option 3: largest generator
            [~,idx] = max(mpc.gen(:,PMAX).*mpc.gen(:,GEN_STATUS));
            ref = full(nmap(mpc.gen(idx,1)));
    end
    bidx.ref(ref) = true;
    bidx.pv(ref)  = false;
end