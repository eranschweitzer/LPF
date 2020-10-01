clear varialbes;
close all;
opt = 1;
%%
conn = dbconn();
cases = fetch(conn,'select * from linpf.cases');
caseid = cases.id;
cases = cases.name;
conn.close();
%%
define_constants;
mpopt = mpoption('out.all', 0, 'verbose', 0);
%%
fprintf('----------------------------------------\n')
fprintf('case                  v_rms    pf_rms   qf_rms   pt_rms   qf_rms\n')
fprintf('--------------------  -------  -------  -------  -------  -------\n')
for k = 1:6%caseid.'
    mpc = loadcase(cases{k});
    r   = runpf(mpc,mpopt);
    switch opt
      case 1
        vars = zhigang(mpc);
      case 2
        vars = SLPF(mpc);
    end
    Cv  = eval_criteria(vars.v,r.bus(:,VM));
    switch opt
      case 1
        Cpf = eval_criteria(real(vars.sf), r.branch(:,PF)/r.baseMVA);
        Cpt = eval_criteria(real(vars.st), r.branch(:,PT)/r.baseMVA);
        Cqf = eval_criteria(imag(vars.sf), r.branch(:,QF)/r.baseMVA);
        Cqt = eval_criteria(imag(vars.st), r.branch(:,QT)/r.baseMVA);
      case 2
        Cpf = eval_criteria(vars.p.f, r.branch(:,PF)/r.baseMVA);
        Cpt = eval_criteria(vars.p.t, r.branch(:,PT)/r.baseMVA);
        Cqf = eval_criteria(vars.q.f, r.branch(:,QF)/r.baseMVA);
        Cqt = eval_criteria(vars.q.t, r.branch(:,QT)/r.baseMVA);
    end
    fprintf('%20s  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f\n',cases{k}, Cv.rms, Cpf.rms, Cqf.rms, Cpt.rms, Cqt.rms)
end