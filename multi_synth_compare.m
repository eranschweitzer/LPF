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
fname = {'RT3000_rndslv_Qlims13-02-2018_1249inputstamp_13-02-2018_1219_ind%d.mat'};
fname = {'ER_deg47_Qlims14-02-2018_0834inputstamp_13-02-2018_1801_ind%d.mat'};
fname = {'RT2500_Qlims14-02-2018_1848inputstamp_14-02-2018_1432_ind%d.mat'};
fname = {'polish2383wp_rndslv_pminfix219-02-2018_2112inputstamp_19-02-2018_2030_ind%d.mat'};
% fname = {'random_solve29-01-2018_1527inputstamp_29-01-2018_1525_ind%d.mat',...
%          'polish2383wp_eatest_ind%d.mat'};
% fname = {'polish_rndslv_Qflowlimit30-01-2018_2117inputstamp_30-01-2018_2041_ind%d.mat'};
% fname = {'polish_ea_Qlims11-02-2018_0156inputstamp_10-02-2018_1640_ind%d.mat'};
% fname = {'case118_rndslv_Qflowlimit31-01-2018_1009inputstamp_31-01-2018_1008_ind%d.mat'};
% fname = {'case118_rates.mat'};
% fname = {'case118_Qflowlimit31-01-2018_1049inputstamp_31-01-2018_1048_ind%d.mat'};
% fname = {'case118_Qflowlimit_loss331-01-2018_1100inputstamp_31-01-2018_1056_ind%d.mat'};
% fname = {'case118_NoQflowlimit_loss331-01-2018_1131inputstamp_31-01-2018_1127_ind%d.mat'};
% fname = {'case118_full_Qlim_loss331-01-2018_1229inputstamp_31-01-2018_1228.mat'};
% fname = {'case118_full_loss3_loadcorrect12-02-2018_1241inputstamp_12-02-2018_1240.mat'}; %no good!
ind = 7;
% deg  = [3.5, 3, 2.5];
deg  = [2.3, 2.5, 3, 3.5];
cmp_case = 'case2383wp';
% fname = fname(4);
%% run comparison
% Thist = cell(7,1);
Thist = cell(ind*length(fname),1);
fitness = cell(length(Thist),1);
for f = 1:length(fname)
    for k = 0:ind-1
        mpc = loadcase(fullfile(fpath,sprintf(fname{f},k)));
    %     synth_compare(mpc,'case2383wp','plot',true,'title', sprintf('Ind%d',k));
%         Thist{k+1} = synth_compare(mpc,cmp_case);
        [Thist{(f-1)*ind + k+1}, fitness{(f-1)*ind + k+1}] = synth_compare(mpc,cmp_case);
    end
%     synth_compare_plots(Thist,'title', sprintf('Avg. Deg: %0.1f',deg(f)))
end
synth_compare_plots(Thist,'FontSize',26)

if length(fitness) > 1
    fitness_range = struct();
    for f = fieldnames(fitness{1}.kl).'
        fitness_range.(f{1}) = struct('kl', [inf, 0, 0], 'bc', [inf, 0, 0]);
    end
    for k = 1:length(fitness)
        for f = fieldnames(fitness{k}.kl).'
            for ff = {'kl', 'bc'}
                % min
                fitness_range.(f{1}).(ff{1})(1) = ...
                    min(fitness_range.(f{1}).(ff{1})(1), fitness{k}.(ff{1}).(f{1}));
                %avg
                fitness_range.(f{1}).(ff{1})(2) = ...
                    fitness_range.(f{1}).(ff{1})(2) + fitness{k}.(ff{1}).(f{1})/length(fitness);
                % max
                fitness_range.(f{1}).(ff{1})(3) = ...
                    max(fitness_range.(f{1}).(ff{1})(3), fitness{k}.(ff{1}).(f{1}));
            end
        end
    end 
    % make tex friendly rows
    texrows = [];
    for f = {'delta', 'pf', 'qf', 'sf'}
        texrows = [texrows, sprintf('%0.3f & %0.3f & %0.3f & %0.3f & %0.3f & %0.3f\\\\\n',...
            fitness_range.(f{1}).kl, fitness_range.(f{1}).bc)];
    end
end