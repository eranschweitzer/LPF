function out = avg_fitness_plots(S, idx)
inputs = {...
        struct('prop','delta', 'label', 'Angel Difference [degrees]'),...
        struct('prop','pf',    'label', 'Real Power (from-end) [MW]'),...
        struct('prop', 'qf',   'label', 'Reactive Power (from-end) [MVAr]'),...
        struct('prop', 'sf',   'label', 'Apparent Power (from-end) [MVA]'),...
        struct('prop', 'sil',   'label', 'Fraction of SIL')};
if nargin == 1
    idx = 2;
end

for k = inputs
    data = zeros(length(S),1);
    for gp = 1:length(S)
        data(gp) = S(gp).(k{1}.prop).bc(idx);
    end
    out.(k{1}.prop) = data;
    figure;
    bar(data)
    ylabel('Avg D_H')
    title(k{1}.label)
end