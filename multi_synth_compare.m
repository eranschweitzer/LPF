close all; clear variables;
%% inputs
fpath = '/Users/eran/Dropbox/ASU/SINE/python/parameter_assignment/mpc';
fname = {'ER_rndslv_deg3529-01-2018_2005inputstamp_29-01-2018_1833_ind%d.mat',...
         'ER_rndslv_deg329-01-2018_1821inputstamp_29-01-2018_1730_ind%d.mat',...
         'ER_test27-01-2018_1547inputstamp_26-01-2018_2144_ind%d.mat'};
fname = {'ER_rndslv_deg23_Qlims12-02-2018_1933inputstamp_12-02-2018_1746_ind%d.mat',...
         'ER_rndslv_deg25_Qlims12-02-2018_2157inputstamp_12-02-2018_1741_ind%d.mat',...
         'ER_rndslv_deg3_Qlims12-02-2018_1919inputstamp_12-02-2018_1739_ind%d.mat',...
         'ER_rndslv_deg35_Qlims12-02-2018_2139inputstamp_12-02-2018_1736_ind%d.mat'};
fname = {'ER243_silconstr77_cbalter07-03-2018_0011inputstamp_06-03-2018_1227_ind%d',...
         'ER286_silconstr77_cbalter07-03-2018_1123inputstamp_07-03-2018_0133_ind%d',...
         'ER320_silconstr77_cbalter09-03-2018_0431inputstamp_09-03-2018_0012_ind%d'};
% fname = {'RT2383_silconstr77_cbalter07-03-2018_0938inputstamp_06-03-2018_2225_ind%d',...
%          'ER286_silconstr77_cbalter07-03-2018_1123inputstamp_07-03-2018_0133_ind%d'};
% fname = {'RT3000_rndslv_Qlims13-02-2018_1249inputstamp_13-02-2018_1219_ind%d.mat'};
% fname = {'ER_deg47_Qlims14-02-2018_0834inputstamp_13-02-2018_1801_ind%d.mat'};
% fname = {'RT2500_Qlims14-02-2018_1848inputstamp_14-02-2018_1432_ind%d.mat'};
% fname = {'polish2383wp_rndslv_pminfix219-02-2018_2112inputstamp_19-02-2018_2030_ind%d.mat'};
% fname = {'random_solve29-01-2018_1527inputstamp_29-01-2018_1525_ind%d.mat',...
%          'polish2383wp_eatest_ind%d.mat'};
% fname = {'polish_rndslv_Qflowlimit30-01-2018_2117inputstamp_30-01-2018_2041_ind%d.mat'};
% fname = {'polish_ea_Qlims11-02-2018_0156inputstamp_10-02-2018_1640_ind%d.mat'};
% fname = {'polish2383wp_silconstr_cbalter05-03-2018_0337inputstamp_04-03-2018_1902_ind%d.mat'};
% fname = {'polish2383wp_silconstr77_cbalter06-03-2018_0421inputstamp_05-03-2018_2120_ind%d.mat'};
% fname = {'case118_rndslv_Qflowlimit31-01-2018_1009inputstamp_31-01-2018_1008_ind%d.mat'};
% fname = {'case118_rates.mat'};
% fname = {'case118_Qflowlimit31-01-2018_1049inputstamp_31-01-2018_1048_ind%d.mat'};
% fname = {'case118_Qflowlimit_loss331-01-2018_1100inputstamp_31-01-2018_1056_ind%d.mat'};
% fname = {'case118_NoQflowlimit_loss331-01-2018_1131inputstamp_31-01-2018_1127_ind%d.mat'};
% fname = {'case118_full_Qlim_loss331-01-2018_1229inputstamp_31-01-2018_1228.mat'}; %used in paper
% fname = {'case118_full_loss3_loadcorrect12-02-2018_1241inputstamp_12-02-2018_1240.mat'}; %no good!

%%%%% everything
fname  = {'polish2383wp_silconstr77_cbalter06-03-2018_0421inputstamp_05-03-2018_2120_ind%d.mat',...
          'RT2383_deg260to243_silconstr77_cbalter13-03-2018_0938inputstamp_12-03-2018_1925_ind%d.mat',...
          'RT2383_deg260_silconstr77_cbalter12-03-2018_2305inputstamp_12-03-2018_0955_ind%d.mat',...
          'RT2383_silconstr77_cbalter07-03-2018_0938inputstamp_06-03-2018_2225_ind%d',...
          'ER243_silconstr77_cbalter07-03-2018_0011inputstamp_06-03-2018_1227_ind%d',...
          'ER286_silconstr77_cbalter07-03-2018_1123inputstamp_07-03-2018_0133_ind%d',...
          'ER320_silconstr77_cbalter09-03-2018_0431inputstamp_09-03-2018_0012_ind%d'};
ind = 7;
% deg  = [3.5, 3, 2.5];
deg  = [2.3, 2.5, 3, 3.5];
cmp_case = 'case2383wp';
% fname = fname(1:3);
% s = struct();
% for prop = {'VMAX', 'RATE_A', 'PMAX', 'QMIN', 'QMAX'}
%     s.(prop{1}).type = 'unbnd';
% end
% 
% for prop = {'ANGMIN', 'ANGMAX'}
%     prop = prop{1};
%     s.(prop).type = 'none';
% end
% 
% [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%     VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% mpopt = mpoption('opf.dc.solver', 'MIPS', 'opf.ac.solver', 'MIPS', 'mips.step_control', 1);
%% run comparison
% Thist = cell(7,1);
Thist = cell(ind*length(fname),1);
fitness = cell(length(Thist),1);
keep    = true(length(Thist),1);
groups  = zeros(length(fname),1);
for f = 1:length(fname)
    fprintf('f%d: ', f)
    for k = 0:ind-1
        fprintf('%d ',k)
        mpc = loadcase(fullfile(fpath,sprintf(fname{f},k)));
%         mpc.bus(:,VMIN) = 0.9; mpc.bus(:,VMAX) = 1.1;
% %         s.RATE_A.idx  = (1:size(mpc.branch,1))';
% %         s.RATE_A.ub   = mpc.branch(:,RATE_A);
%         mpc.softlims = s;
%         mpc = toggle_softlims(mpc, 'on');
%         r = runopf(mpc, mpopt);
    %     synth_compare(mpc,'case2383wp','plot',true,'title', sprintf('Ind%d',k));
%         Thist{k+1} = synth_compare(mpc,cmp_case);
        idx = (f-1)*ind + k+1;
        if idx == 1
            [Thist{idx}, fitness{idx}, ~, ~, mpcref] = synth_compare(mpc,cmp_case);
        else
            [Thist{idx}, fitness{idx}] = synth_compare(mpc,mpcref);
        end
        if ~isstruct(Thist{idx})
            warning('removing individual %d (file %s)', k, sprintf(fname{f}, k))
            keep(idx) = false;
        else
            groups(f) = groups(f) + 1;
        end
    end
    fprintf('\n')
%     synth_compare_plots(Thist,'title', sprintf('Avg. Deg: %0.1f',deg(f)))
end
Thist = Thist(keep);
fitness = fitness(keep);
% return
synth_compare_plots(Thist,'FontSize',26, 'groups', groups)

if length(fitness) > 1
    groupcnt = [0; cumsum(groups)];
    fitness_range(length(groups)) = struct();
    for gp = 1:length(groups) 
        for f = fieldnames(fitness{1}.kl).'
            fitness_range(gp).(f{1}) = struct('kl', [inf, 0, 0], 'bc', [inf, 0, 0]);
        end
        for k = groupcnt(gp)+1:groupcnt(gp+1)%1:length(fitness)
            for f = fieldnames(fitness{k}.kl).'
                for ff = {'kl', 'bc'}
                    % min
                    fitness_range(gp).(f{1}).(ff{1})(1) = ...
                        min(fitness_range(gp).(f{1}).(ff{1})(1), fitness{k}.(ff{1}).(f{1}));
                    %avg
                    fitness_range(gp).(f{1}).(ff{1})(2) = ...
                        fitness_range(gp).(f{1}).(ff{1})(2) + fitness{k}.(ff{1}).(f{1})/groups(gp);%length(fitness);
                    % max
                    fitness_range(gp).(f{1}).(ff{1})(3) = ...
                        max(fitness_range(gp).(f{1}).(ff{1})(3), fitness{k}.(ff{1}).(f{1}));
                end
            end
        end 
        texrows = [];
        for f = {'delta', 'pf', 'qf', 'sf', 'sil'}
            texrows = [texrows, sprintf('%0.3f & %0.3f & %0.3f\\\\\n', fitness_range(gp).(f{1}).bc)];
        end
        fitness_range(gp).texrows = texrows;
    end
%     % make tex friendly rows
%     texrows = [];
%     for f = {'delta', 'pf', 'qf', 'sf', 'sil'}
%         texrows = [texrows, sprintf('%0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f\\\\\n',...
%             fitness_range.(f{1}).kl, fitness_range.(f{1}).bc)];
%     end
end

%% save csv
% for k = 1:ind
%     writetable(mystruct2table(Thist{k}),sprintf('~/Dropbox/Publications/2018parameter_assignment/figs/ER47_%d.csv',k));
% end