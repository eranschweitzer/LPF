clear variables; close all;
%% load case
tablename = 'linpf.res';
batchsize = 50;
precheck = true;
define_constants;
casename = 'case118';
mpc = loadcase(casename);

N = size(mpc.bus,1);
M = size(mpc.branch,1);
G = size(mpc.gen,1);
baseMVA = mpc.baseMVA;
%% maps and operating matrices
nmap  = sparse(mpc.bus(:,BUS_I),1,1:N);
E  = sparse([1:M,1:M].',[full(nmap(mpc.branch(:,F_BUS)));full(nmap(mpc.branch(:,T_BUS)))],[ones(M,1);-1*ones(M,1)],M,N);
F  = sparse(1:M,full(nmap(mpc.branch(:,F_BUS))),1,M,N);
T  = sparse(1:M,full(nmap(mpc.branch(:,T_BUS))),1,M,N);
gbus  = full(nmap(mpc.gen(:,GEN_BUS))); %generator buses
gmap  = sparse(gbus,1:G,(mpc.gen(:,GEN_STATUS) > 0),N,G); %maps generators onto buses
bus_with_ongen = sum(gmap,2) > 0;

%% bus types
[pq_idx,pv_idx,ref_idx] = bustype_map(mpc,bus_with_ongen);
Npq = sum(pq_idx); Npv = sum(pv_idx); Nref= sum(ref_idx);

%% u0
vg = ones(N,1); vg(gbus) = mpc.gen(:,VG); 
u0 = u0init(vg,pv_idx,ref_idx);
theta_ref = mpc.bus(ref_idx,VA)*pi/180;

%% matrices fixing PV and Ref quantities
[matI,matU] = pv_ref_mats(pv_idx,ref_idx,Npv,Nref,N,u0);
%% load and generation
Sp = powervectors(mpc,gmap);

%% branch parts
Sb = BranchParts(mpc);

%% AC solution
mpcac  = runpf(mpc);
vtrue  = mpcac.bus(:,VM);        %true voltage magnitude solution
ttrue  = mpcac.bus(:,VA)*pi/180; %true voltage angle solution
Pftrue = mpcac.branch(:,PF);     %true real power at from bus [MW]
Pttrue = mpcac.branch(:,PT);     %true real power at to bus [MW]
Qftrue = mpcac.branch(:,QF);     %true reactive power at from bus [MVAr]
Qttrue = mpcac.branch(:,QT);     %true reactive power at to bus [MVAr]
%% iteration IDs
IDs = genIDs();
Ntotal = size(IDs,1);
colnames = {'casename', 'id', 'prop', 'criteria', 'value'};
missedID = int8(zeros(0, size(IDs,2)));
%% Main Loop

parpool(60);
parfor batch = 1:ceil(Ntotal/batchsize)
    IDbatch = IDs(batchsize*(batch-1)+1:min(batchsize*batch,Ntotal),:);
    if precheck
        if checkbatch(casename,IDbatch)
            continue
        end
    end
    Tres = cell(size(IDbatch,1),1);
    for k = 1:size(IDbatch,1)
        ids  = unpackIDs(IDbatch(k,:));
        vars = pfsolve(ids,F,T,E,Sb,Sp,matI,matU,theta_ref,Nref,Npv,N,u0);

        % Evaluation criteria for voltage and angle 
        Cv   = eval_criteria(vars.v,vtrue);
        Ct   = eval_criteria(vars.theta,ttrue);

        % Evalutaion criteria for flow calculation
        Cpf = cell(6,2); Cpt = cell(6,2);
        Cqf = cell(6,2); Cqt = cell(6,2);
        for ID = 0:5
            for btype = 0:1
                P = calcPflow(ID,btype,vars,F,T,E,Sb);
                Q = calcQflow(ID,btype,vars,F,T,E,Sb);
                Cpf{ID+1,btype+1} = eval_criteria(baseMVA*P.f,Pftrue);
                Cpt{ID+1,btype+1} = eval_criteria(baseMVA*P.t,Pttrue);
                Cqf{ID+1,btype+1} = eval_criteria(baseMVA*Q.f,Qftrue);
                Cqt{ID+1,btype+1} = eval_criteria(baseMVA*Q.t,Qttrue);
            end
        end
        
        % create result table 
        Tres{k} = result2table(colnames,casename,ids,vars.convg,Cv,Ct,Cpf,Cpt,Cqf,Cqt)
    end
    % write to database
    flag = result2db(tablename,colnames,Tres)
    
    % if write to database didn't work save the the missing IDs for later
    % handling
    if ~flag
        missedID = [missedID;IDbatch];
    end
end
delete(gcp);

%% save missed IDs
if size(missedID,1) > 0
    save([casename, '_missedIDs.mat'], 'missedID')
end