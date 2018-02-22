function Cf = flow_performance(flids, ids, vars, vtrue, F, T, E, Sb, baseMVA)

Cpf = cell(length(flids.p),1); Cpt = cell(length(flids.p),1);
Cqf = cell(length(flids.q),1); Cqt = cell(length(flids.q),1);
for k = 1:length(flids.p)
    flid = str2ids(flids.p{k});
    flid.c = ids.c;
    if flid.Q == 3
    end
    P = flowcalc(flid, 'real', vars, F, T, E, Sb);
    Cpf{k} = eval_criteria(P.f, vtrue.pf/baseMVA);
    Cpt{k} = eval_criteria(P.t, vtrue.pt/baseMVA);
end
for k = 1:length(flids.q)
    flid = str2ids(flids.q{k});
    flid.c = ids.c;
    if flid.Q == 3
    end
    Q = flowcalc(flid, 'imag', vars, F, T, E, Sb);
    Cqf{k} = eval_criteria(Q.f, vtrue.qf/baseMVA);
    Cqt{k} = eval_criteria(Q.t, vtrue.qt/baseMVA);
end

Cf = struct('pf', {Cpf}, 'qf', {Cqf}, 'pt', {Cpt}, 'qt', {Cqt});