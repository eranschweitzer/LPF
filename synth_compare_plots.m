function synth_compare_plots(Thist,varargin)
plt_title  = varargin_parse(varargin,'title', '');
fnt_size   = varargin_parse(varargin,'FontSize', 16);
groupcnt  = varargin_parse(varargin,'groups',0);

cmap = [0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840
        0         0.4470    0.7410];
        
ind = numel(Thist);
if (ind == 1) && iscell(Thist)
    Thist = Thist{1};
end
if groupcnt ~= 0
   groupcnt = cumsum(groupcnt);
else
   groupcnt = ind;
end

inputs = {...
        struct('prop','delta', 'label', 'Angel Difference [degrees]'),...
        struct('prop','pf',    'label', 'Real Power (from-end) [MW]'),...
        struct('prop', 'qf',   'label', 'Reactive Power (from-end) [MVAr]'),...
        struct('prop', 'sf',   'label', 'Apparent Power (from-end) [MVA]'),...
        struct('prop', 'sil',   'label', 'Fraction of SIL')};
    
for k = inputs
    figure;
    subplot(1,2,1)
    xl     = sprintf('%scenters',k{1}.prop);
    ylmpc  = sprintf('%smpc', k{1}.prop);
    ylsynth= sprintf('%ssynth', k{1}.prop);
    if ind >  1
        plot(Thist{1}.(xl),Thist{1}.(ylmpc),'-ok','linewidth', 3);
    else
        plot(Thist.(xl),Thist.(ylmpc),'-o','linewidth', 3);
    end
    hold on;
    if ind >  1
        for i = 1:ind
            c = find(groupcnt >= i, 1);
            plot(Thist{i}.(xl),Thist{i}.(ylsynth),'-*','color',cmap(c,:));
        end
    else
        plot(Thist.(xl),Thist.(ylsynth),'-*');
    end
    set(gca,'Yscale','log')
    set(gca,'FontSize',fnt_size);
    xlabel(k{1}.label)
    ylabel('Density')
    if ~strcmp(plt_title,'')
        title(plt_title)
    end

    subplot(1,2,2)
    ylmpc  = sprintf('%smpcq', k{1}.prop);
    ylsynth= sprintf('%ssynthq', k{1}.prop);
    ymax = 0;
    if ind >  1
        plot(Thist{1}.q,Thist{1}.(ylmpc),'k','linewidth', 3);
        ymax = max(ymax,max(Thist{1}.(ylmpc)));
    else
        plot(Thist.q,Thist.(ylmpc),'linewidth', 3);
        ymax = max(ymax,max(Thist.(ylmpc)));
    end
    hold on;
    if ind >  1
        for i = 1:ind
            c = find(groupcnt >= i, 1);
            plot(Thist{i}.q,Thist{i}.(ylsynth),'--','color',cmap(c,:));
            ymax = max(ymax,max(Thist{i}.(ylsynth)));
        end
    else
        plot(Thist.q,Thist.(ylsynth));
        ymax = max(ymax,max(Thist.(ylsynth)));
    end
    set(gca,'Yscale','log')
    set(gca,'FontSize',fnt_size);
    ylabel(k{1}.label)
    xlabel('quantile')
    if ~strcmp(plt_title,'')
        title(plt_title)
    end
    ax = gca;
    p = 0;
    while ymax/10^p > 1
        p = p + 1;
    end
    ax.YLim(2) = 10^p;
    ax.YLim(1) = 1e-5;
end

if length(groupcnt) > 1
    figure;
    title('legend')
    for k = 1:length(groupcnt)
        plot([0 1], [k k], 'color', cmap(k,:))
        text(0.5,k, sprintf('group %d', k))
        if k == 1
            hold on;
        end
    end
end
%% Angle
% figure;
% subplot(1,2,1)
% if ind >  1
%     plot(Thist{1}.deltacenters,Thist{1}.deltampc,'-o','linewidth', 3);
% else
%     plot(Thist.deltacenters,Thist.deltampc,'-o','linewidth', 3);
% end
% hold on;
% if ind >  1
%     for k = 1:ind
%         plot(Thist{k}.deltacenters,Thist{k}.deltasynth,'-*k');
%     end
% else
%     plot(Thist.deltacenters,Thist.deltasynth,'-*');
% end
% set(gca,'Yscale','log')
% xlabel('Angel Difference [degrees]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% if ind >  1
%     plot(Thist{1}.q,Thist{1}.deltampcq,'linewidth', 3);
% else
%     plot(Thist.q,Thist.deltampcq,'linewidth', 3);
% end
% hold on;
% if ind >  1
%     for k = 1:ind
%         plot(Thist{k}.q,Thist{k}.deltasynthq);
%     end
% else
%     plot(Thist.q,Thist.deltasynthq);
% end
% set(gca,'Yscale','log')
% ylabel('Angel Difference [degrees]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;
% 
% %% Real Powerflow (from end)
% figure;
% subplot(1,2,1)
% plot(Thist.pfcenters,Thist.pfmpc,'-o','linewidth', 3);
% hold on;
% plot(Thist.pfcenters,Thist.pfsynth,'-*');
% set(gca,'Yscale','log')
% xlabel('Real Power (from-end) [MW]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% plot(Thist.q,Thist.pfmpcq,'linewidth', 3);
% hold on;
% plot(Thist.q,Thist.pfsynthq);
% set(gca,'Yscale','log')
% ylabel('Real Power (from-end) [MW]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;
% 
% %% Reactive Powerflow (from end)
% figure;
% subplot(1,2,1)
% plot(Thist.qfcenters,Thist.qfmpc,'-o','linewidth', 3);
% hold on;
% plot(Thist.qfcenters,Thist.qfsynth,'-*');
% set(gca,'Yscale','log')
% xlabel('Reactive Power (from-end) [MVAr]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% plot(Thist.q,Thist.qfmpcq,'linewidth', 3);
% hold on;
% plot(Thist.q,Thist.qfsynthq);
% set(gca,'Yscale','log')
% ylabel('Reactive Power (from-end) [MVAr]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;
% 
% %% Apparent Powerflow (from end)
% figure;
% subplot(1,2,1)
% plot(Thist.sfcenters,Thist.sfmpc,'-o','linewidth', 3);
% hold on;
% plot(Thist.sfcenters,Thist.sfsynth,'-*');
% set(gca,'Yscale','log')
% xlabel('Apparent Power (from-end) [MVA]')
% ylabel('Density')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% 
% subplot(1,2,2)
% plot(Thist.q,Thist.sfmpcq,'linewidth', 3);
% hold on;
% plot(Thist.q,Thist.sfsynthq);
% set(gca,'Yscale','log')
% ylabel('Apparent Power (from-end) [MVA]')
% xlabel('quantile')
% if ~strcmp(plt_title,'')
%     title(plt_title)
% end
% ax = gca;
% ax.YLim(1) = 1e-5;