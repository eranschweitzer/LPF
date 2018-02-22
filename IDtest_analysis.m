%% Analysis of linpf test to determine optimal ID
clear variables; close all;
%% Database connection
conn = dbconn();
best_table = 'linpf.best0';
table_name = 'linpf.res';
actest_table = 'linpf.actest';
iterID = 0; %this should match the best_table choice
plots = true;
%%
casenames = fetch(conn,['SELECT distinct casename FROM ', best_table, ' ORDER BY casename']);
casenames = casenames.casename;
%% Determine v* (optimal values)
vstar = struct('volt', zeros(length(casenames),3), 'ang', zeros(length(casenames),3));
cols = {'cor', 'max', 'avg'};
for x = {'volt', 'ang'}
    for c = {'max', 'avg'}
        sql = ['SELECT casename, min(value) as v FROM ', best_table,...
                ' WHERE prop=''', x{1},''' and criteria=''',c{1},''' GROUP BY casename ORDER BY casename'];
        data = fetch(conn,sql);
        vstar.(x{1})(:,strcmp(cols,c{1})) = data.v;
    end
    sql = ['SELECT casename, max(value) as v FROM ', best_table,...
                ' WHERE prop=''', x{1},''' and criteria=''cor'' GROUP BY casename ORDER BY casename'];
    data = fetch(conn,sql);
    vstar.(x{1})(:,strcmp(cols,'cor')) = data.v;
end

vstar.volt = array2table(vstar.volt,'VariableNames',cols,'RowNames',casenames);
vstar.ang = array2table(vstar.ang,'VariableNames',cols,'RowNames',casenames);
%%
IDs = genIDs('iterID',iterID); 
res = struct('id',cell(1,1), 'cost', []);
for d = 1:size(IDs,1)
    if mod(d,1000) == 0
        fprintf('Working on ID %d of %d\n', d, size(IDs,1))
    end
    idstr = ids2str(unpackIDs(IDs(d,:)));
    sql = ['SELECT * FROM ', table_name,...
           ' WHERE id=''', idstr, ''' and prop in (''volt'', ''ang'') ',...
           'and criteria in (''cor'', ''max'', ''avg'') ',...
           'ORDER BY casename, prop, criteria'];
    data = fetch(conn,sql);
    if length(data.value) ~= 2*3*length(casenames)
        % this ID did not converge for all cases, do not include
        continue
    end
    % test that using this id all cases converged for the AC powerflow
    sql = ['SELECT linpf FROM ', actest_table, ' WHERE id=''', idstr, ''''];
    actest = fetch(conn,sql);
    if ~all(actest.linpf)
        continue
    end
    cost = 0;
    for i = casenames.'
        for x = {'volt', 'ang'}
            for c = cols
                mask = strcmp(data.casename,i{1}) & ...
                       strcmp(data.prop, x{1})    & ...
                       strcmp(data.criteria, c{1});
                cost = cost + (vstar.(x{1}).(c{1})(i{1}) - data.value(mask)).^2;
            end
        end
    end
    res.id{end+1} = idstr;
    res.cost(end+1) = cost;
end

%% get minimum cost
[~,dstar] = min(res.cost);

%% get data
sql = ['SELECT * FROM ', table_name,...
           ' WHERE id=''', res.id{dstar}, ''' and prop in (''volt'', ''ang'') ',...
           'and criteria in (''cor'', ''max'', ''avg'') ',...
           'ORDER BY casename, prop, criteria'];
data = fetch(conn,sql);

%% output
out = vstar;
for x = {'volt', 'ang'}
    tmp = out.(x{1});
    for i = casenames.'
        for c = cols
            mask = strcmp(data.casename,i{1}) & ...
                       strcmp(data.prop, x{1})    & ...
                       strcmp(data.criteria, c{1});
            tmp.(c{1})(i{1}) = data.value(mask);
        end 
    end
    tmp.Properties.VariableNames = {'corr_actual', 'max_actual', 'avg_actual'};
    out.(x{1}) = [out.(x{1}) tmp];
end
%% plots
if ~plots
    return
end
figure;
subplot(1,3,1)
bar(categorical(out.volt.Properties.RowNames), [out.volt.cor, out.volt.corr_actual])
ylabel('v correlation')

subplot(1,3,2)
bar(categorical(out.volt.Properties.RowNames), [out.volt.max, out.volt.max_actual])
ylabel('v max')

subplot(1,3,3)
bar(categorical(out.volt.Properties.RowNames), [out.volt.avg, out.volt.avg_actual])
ylabel('v avg')

figure;
subplot(1,3,1)
bar(categorical(out.ang.Properties.RowNames), [out.ang.cor, out.ang.corr_actual])
ylabel('\theta correlation')

subplot(1,3,2)
bar(categorical(out.ang.Properties.RowNames), [out.ang.max, out.ang.max_actual])
ylabel('\theta max')
set(gca, 'YScale', 'log')

subplot(1,3,3)
bar(categorical(out.ang.Properties.RowNames), [out.ang.avg, out.ang.avg_actual])
ylabel('\theta avg')
set(gca, 'YScale', 'log')