function PFbenchmark(casename,varargin)

%% optional inputs
precheck = varargin_parse(varargin,'precheck',true);
batchsize= varargin_parse(varargin,'batchsize',50);
itermax  = varargin_parse(varargin,'itermax',50);
GwAru    = varargin_parse(varargin,'GwAru','default');
GwAmu    = varargin_parse(varargin,'GWAmu','default');
bcmode   = varargin_parse(varargin,'bcmode','default');
bshmode  = varargin_parse(varargin,'bshmode','default');
udefault = varargin_parse(varargin,'udefault','uref');
uopt     = varargin_parse(varargin,'uopt','none');
tablename= varargin_parse(varargin,'tablename','linpf.res');
testtype = varargin_parse(varargin,'testtype','default');
caseid   = varargin_parse(varargin,'caseid', 0);
%%
fprintf('---------------------------\n')
fprintf('Running case: %s\n',casename)
fprintf('---------------------------\n')
%% load case
define_constants;
mpc = loadcase(casename);

N = struct('t',size(mpc.bus,1));
M = size(mpc.branch,1);
G = size(mpc.gen,1);
baseMVA = mpc.baseMVA;
%% maps and operating matrices
nmap  = sparse(mpc.bus(:,BUS_I),1,1:N.t);
E  = sparse([1:M,1:M].',[full(nmap(mpc.branch(:,F_BUS)));full(nmap(mpc.branch(:,T_BUS)))],[ones(M,1);-1*ones(M,1)],M,N.t);
F  = sparse(1:M,full(nmap(mpc.branch(:,F_BUS))),1,M,N.t);
T  = sparse(1:M,full(nmap(mpc.branch(:,T_BUS))),1,M,N.t);
gbus  = full(nmap(mpc.gen(:,GEN_BUS))); %generator buses
gstat = mpc.gen(:,GEN_STATUS) > 0;  % mask of dispatched generators
gmap  = sparse(gbus,1:G,gstat,N.t,G); %maps generators onto buses
bus_with_ongen = sum(gmap,2) > 0;

%% bus types
bidx = bustype_map(mpc,bus_with_ongen,gstat,nmap);
N.pq = sum(bidx.pq); N.pv = sum(bidx.pv); N.ref= sum(bidx.ref);

%% u0
vg = ones(N.t,1); vg(gbus(gstat)) = mpc.gen(gstat,VG); 
u0 = u0init(vg,bidx,'udefault',udefault);
theta_ref = mpc.bus(bidx.ref,VA)*pi/180;

%% matrices fixing PV and Ref quantities
% [matI,matU] = pv_ref_mats(pv_idx,ref_idx,Npv,Nref,N,u0);
%% load and generation
Sp = powervectors(mpc,gmap);
Sp.Qg(bidx.pv) = 0;
Sp.Qg(bidx.ref) = 0;
Sp.Pg(bidx.ref) = 0;
%% branch parts
Sb = BranchParts(mpc);
Gw = branchweights(Sb,'GwAru',GwAru,'GwAmu',GwAmu,'bcmode',bcmode,'bshmode',bshmode);
Sb.b    = Gw.b*Sb.b;
Sb.balt = Gw.b*Sb.balt;
Sb.bc   = Gw.bc*Sb.bc;
Sb.bsh  = Gw.bsh*Sb.bsh;
Sb.g    = Gw.g*Sb.g;
%% AC solution
mpopt   = mpoption; mpopt.out.all = 0;
mpcac  = runpf(mpc,mpopt);
vtrue    = struct();
vtrue.v  = mpcac.bus(:,VM);        %true voltage magnitude solution
vtrue.t  = mpcac.bus(:,VA);        %true voltage angle solution [Deg]
vtrue.pf = mpcac.branch(:,PF);     %true real power at from bus [MW]
vtrue.pt = mpcac.branch(:,PT);     %true real power at to bus [MW]
vtrue.qf = mpcac.branch(:,QF);     %true reactive power at from bus [MVAr]
vtrue.qt = mpcac.branch(:,QT);     %true reactive power at to bus [MVAr]
%% iteration IDs
if strcmp(testtype,'zhcomp')
    IDs = newIDs();
    colnames = struct('vt', {{'caseid', 'id', 'prop', 'criteria', 'value'}},...
                      'flow', {{'caseid', 'pfid', 'id', 'prop', 'criteria', 'value'}});
%     missedID = cell(length(IDs),1);
else
    IDs = genIDs();
    colnames = {'casename', 'id', 'prop', 'criteria', 'value'};
%     missedID = int8(zeros(0, size(IDs,2)));    
end
acconvg_names = {'casename','id','dc','linpf'};
Ntotal = size(IDs,1);


%% Main Loop

% parpool(60);
for batch = 1:ceil(Ntotal/batchsize)
    batchrange = batchsize*(batch-1)+1:min(batchsize*batch,Ntotal);
    fprintf('Batch %d. range %d:%d\n', batch, batchrange(1), batchrange(end))
    IDbatch = IDs(batchrange,:); %#ok<PFBNS>
    if precheck
        if checkbatch(tablename,casename,IDbatch)
            continue
        end
    end
%     if strcmp(testtype,'default') || strcmp(testtype,'acconvg')
        Tres = cell(size(IDbatch,1),1);
        Tres2= cell(size(IDbatch,1),1);
        fail = false(size(IDbatch,1),1);
%     end
    for k = 1:size(IDbatch,1)
        ids  = unpackIDs(IDbatch(k,:));
        %fprintf('%d id: %s\n', k, ids2str(ids))
        vars = pfsolve(ids,F,T,E,Sb,Sp,bidx,theta_ref,N,u0,...
               'itermax',itermax,'uopt',uopt);
        if vars.convg == 1
            if strcmp(testtype,'default')
                % Evaluation criteria for voltage and angle (in degrees)
                Cv   = eval_criteria(vars.v,vtrue.v); %#ok<PFBNS>
                Ct   = eval_criteria(vars.theta*180/pi,vtrue.t);

                % Evalutaion criteria for flow calculation
                pfids = [0,1,3,5];
                Cpf = cell(6,2); Cpt = cell(6,2);
                Cqf = cell(6,2); Cqt = cell(6,2);
                for ID = pfids%0:5
                    for btype = 0:1
                        P = calcPflow(ID,btype,vars,F,T,E,Sb);
                        Q = calcQflow(ID,btype,vars,F,T,E,Sb);
                        Cpf{ID+1,btype+1} = eval_criteria(baseMVA*P.f,vtrue.pf);
                        Cpt{ID+1,btype+1} = eval_criteria(baseMVA*P.t,vtrue.pt);
                        Cqf{ID+1,btype+1} = eval_criteria(baseMVA*Q.f,vtrue.qf);
                        Cqt{ID+1,btype+1} = eval_criteria(baseMVA*Q.t,vtrue.qt);
                    end
                end

                % create result table 
                Tres{k} = result2table(colnames,casename,ids,vars.convg,Cv,Ct,Cpf,Cpt,Cqf,Cqt,pfids);
            elseif strcmp(testtype,'zhcomp')
                if caseid == 0
                    error('caseid cannot be 0!')
                end
                % Evaluation criteria for voltage and angle (in radians)
                Cv   = eval_criteria(vars.v,vtrue.v);
                Ct   = eval_criteria(vars.theta, vtrue.t*pi/180);
                Ctd  = eval_criteria(E*vars.theta, E*vtrue.t*pi/180);
                flids = flowids();
                Cf = flow_performance(flids, ids, vars, vtrue, F, T, E, Sb, baseMVA);
                
%                 Tres = struct();
                Tres{k}    = resprep_vt(ids, Cv, Ct, Ctd, caseid);
                Tres2{k}   = resprep_flow(ids,flids,Cf,caseid);
%                 result2db(tablename, colnames, Tres);
%                 if ~flag
%                     missedID{batchsize*(batch-1)+1 + k} = ids2str(ids);
%                 end
            elseif strcmp(testtype,'acconvg')
                actest = linpf2ac_convergence(mpc,vars,bidx);
                Tres{k} = cell2table({casename,ids2str(ids),logical(actest.dc), logical(actest.linpf)},'VariableNames',acconvg_names);   
            else
                error('Unknown test type')
            end
        else
            fail(k) = true;
        end
    end
    fprintf('batch finished.\n')
%     if strcmp(testtype,'default') || strcmp(testtype,'acconvg')
        Tres = Tres(~fail);
        Tres2= Tres2(~fail);
        if ~isempty(Tres)
            % write to database
            if strcmp(testtype,'default')
                flag = result2db(tablename,colnames,Tres);
            elseif strcmp(testtype,'acconvg')
                flag = result2db(tablename,acconvg_names,Tres);
            elseif strcmp(testtype,'zhcomp')
                result2db(tablename.vt, colnames.vt, Tres);
                result2db(tablename.flow, colnames.flow, Tres2);
            end
        else
            % if all diverged,set flag to true so no missed IDs are stored
            flag = true;
        end

        % if write to database didn't work save the the missing IDs for later
        % handling
%         if ~flag
%             missedID = [missedID;IDbatch];
%         end
%     end
end
% delete(gcp);

%% save missed IDs
% try
%     if strcmp(testtype,'zhcomp')
%         missedID = missedID(~cellfun(@isempty,missedID));
%     end
%     if size(missedID,1) > 0
%         save([casename, '_missedIDs.mat'], 'missedID')
%     end
% catch 
% end