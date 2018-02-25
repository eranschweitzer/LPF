clear variables; close all;
%%
batchsize = 1968/24;
precheck = false;
tablename= struct('vt', 'linpf.vt', 'flow', 'linpf.flows');
GwAru    = 'default';
testtype = 'zhcomp';
% cases = {'case118','case14','case24_ieee_rts','case_ieee30',...
%          'case39','case300','case3375wp','case_ACTIVSg2000','case6470rte',...
%          'case_ACTIVSg10k','case13659pegase'};
conn = dbconn();
cases = fetch(conn,'select * from linpf.cases');
caseid = cases.id;
cases = cases.name;
conn.close();
% parpool(24);
%%
for c = 1:length(cases)
    PFbenchmark(cases{c},'precheck',precheck,'GwAru',GwAru,'tablename',tablename,'testtype',testtype, ...
        'caseid', caseid(c), 'batchsize', batchsize)
end

%%
% delete(gcp);