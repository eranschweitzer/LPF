% diagnostics
bid = 2019; % branch where error occurs
% bid = 214;
levels = 3;
Lap = (F - T)'*(F-T);
degmat = diag(Lap);
Adj = diags(degmat) - Lap;

[tmp, flist] = find(F);
[~,order] = sort(tmp);
flist = flist(order);
[tmp, tlist] = find(T);
[~,order] = sort(tmp);
tlist = tlist(order);
%% branch ends
fend = flist(bid);
tend = tlist(bid);

%% BFS
neighbors = sparse([fend;tend], 1, 1, N.t,1);
for k = 1:levels
    neighbors = neighbors | Adj'*neighbors; %BFS
end

nodeids = find(neighbors);
edgeids = find(ismember(flist, nodeids) & ismember(tlist, nodeids));

pqmask  = bidx.pq(nodeids);
pvmask  = bidx.pv(nodeids);
rfmask  = bidx.ref(nodeids);

%% node and branch parameters
% r =  Sb.g./(Sb.g.^2 + Sb.b.^2);
% x = -Sb.b./(Sb.g.^2 + Sb.b.^2);
delta   = E*vars.theta;
bcharge = F'*Sb.bc/2 + T'*Sb.bc/2;

%% graph of neighborhood
tmpmap  = sparse(nodeids, 1, 1:length(nodeids));
% tmpadj  = sparse(full(tmpmap([flist(edgeids);tlist(edgeids)])), ...
%                  full(tmpmap([tlist(edgeids);flist(edgeids)])), 1);
tmpadj  = sparse(full(tmpmap(flist(edgeids))), full(tmpmap(tlist(edgeids))), 1, length(nodeids), length(nodeids));
Grph = digraph(tmpadj>0,cellstr(num2str(nodeids)));
edgelabels = cell(size(Grph.Edges,1),1);
weights = struct('g', zeros(size(Grph.Edges,1),1), 'b', zeros(size(Grph.Edges,1),1),...
    'tau', ones(size(Grph.Edges,1),1), 'tshift', zeros(size(Grph.Edges,1),1),...
    'deltaerr', zeros(size(Grph.Edges,1),1),...
    'bc', bcharge(nodeids), 'bsh',Sb.bsh(nodeids), 'gsh', Sb.gsh(nodeids),...
    'Pd', Sp.Pd(nodeids), 'Qd', Sp.Qd(nodeids),...
    'res', residual(nodeids), 'verr', vars.v(nodeids) - vtrue.v(nodeids));

for k = 1:length(edgeids)
    tmpid = Grph.findedge(full(tmpmap(flist(edgeids(k)))), full(tmpmap(tlist(edgeids(k)))));
    if isempty(edgelabels{tmpid})
        edgelabels{tmpid} = num2str(edgeids(k));
        weights.g(tmpid) = Sb.g(edgeids(k));
        weights.b(tmpid) = Sb.b(edgeids(k));
        weights.tau(tmpid) = Sb.tau(edgeids(k));
        weights.tshift(tmpid) = Sb.tshift(edgeids(k));
        weights.deltaerr(tmpid) = delta(edgeids(k)) - vtrue.delta(edgeids(k));
    else
        edgelabels{tmpid} = sprintf('%s,%d',edgelabels{tmpid},edgeids(k));
        weights.g(tmpid) = weights.g(tmpid) + Sb.g(edgeids(k));
        weights.b(tmpid) = weights.b(tmpid) + Sb.b(edgeids(k));
        if weights.tau(tmpid) ~= Sb.tau(edgeids(k))
            warning('mismatching tau for parallel lines, using tau1*tau2')
            weights.tau(tmpid) = weights.tau(tmpid)*Sb.tau(edgeids(k));
        end
        if weights.tshift(tmpid) ~= Sb.tshift(edgeids(k))
            error('mismatching tshift for parallel lines')
        end
        if abs(weights.deltaerr(tmpid)) ~= abs(delta(edgeids(k)) - vtrue.delta(edgeids(k)))
            error('Two different angle difference errors found on parallel branches.')
        end
    end
end
weights.r =  weights.g./(weights.g.^2 + weights.b.^2);
weights.x = -weights.b./(weights.g.^2 + weights.b.^2);
weights.rx = weights.r./weights.x;
weights.tau = log(weights.tau);
figure;
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
highlight(h,full(tmpmap(fend)),full(tmpmap(tend)))
highlight(h,full(tmpmap(fend)),full(tmpmap(tend)), 'EdgeColor', 'r')
highlight(h,full(tmpmap(nodeids(pqmask))))
highlight(h,full(tmpmap(nodeids(pvmask))), 'NodeColor', 'g', 'MarkerSize',10)
highlight(h,full(tmpmap(nodeids(rfmask))), 'NodeColor', 'k')

%%
figure;
subplot(5,2,1) 
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('r')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'r');
h.LineWidth = w;
h.EdgeColor = repmat(h.EdgeColor, length(w), 1);
h.EdgeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.EdgeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,3) 
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('x')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'x');
h.LineWidth = w;
h.EdgeColor = repmat(h.EdgeColor, length(w), 1);
h.EdgeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.EdgeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,5) 
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('r/x')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'rx');
h.LineWidth = w;
h.EdgeColor = repmat(h.EdgeColor, length(w), 1);
h.EdgeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.EdgeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,7) 
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('log(\tau)')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'tau');
h.LineWidth = w;
h.EdgeColor = repmat(h.EdgeColor, length(w), 1);
h.EdgeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.EdgeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,9)
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('\delta error')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'deltaerr');
h.LineWidth = w;
h.EdgeColor = repmat(h.EdgeColor, length(w), 1);
h.EdgeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.EdgeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

%%% nodal props
subplot(5,2,2)
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('line charging')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'bc', 'vmax', 12, 'vmin', 5);
h.MarkerSize = w;
h.NodeColor = repmat(h.NodeColor, length(w), 1);
h.NodeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.NodeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,4)
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('Bsh')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'bsh', 'vmax', 12, 'vmin', 5);
h.MarkerSize = w;
h.NodeColor = repmat(h.NodeColor, length(w), 1);
h.NodeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.NodeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,6)
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('P load')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'Pd', 'vmax', 12, 'vmin', 5);
h.MarkerSize = w;
h.NodeColor = repmat(h.NodeColor, length(w), 1);
h.NodeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.NodeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,8)
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('Q load')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'Qd', 'vmax', 12, 'vmin', 5);
h.MarkerSize = w;
h.NodeColor = repmat(h.NodeColor, length(w), 1);
h.NodeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.NodeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);

subplot(5,2,10)
h = plot(Grph,'Layout','layered', 'Sources',[full(tmpmap(fend)),full(tmpmap(tend))], 'edgelabel',edgelabels);
title('Residual Error')
set(gca, 'Xtick',[], 'Ytick',[]);
[w,posmask,negmask] = weights_for_graph(weights,'res', 'vmax', 12, 'vmin', 5);
h.MarkerSize = w;
h.NodeColor = repmat(h.NodeColor, length(w), 1);
h.NodeColor(posmask,:) = repmat([0 0 0], sum(posmask), 1);
h.NodeColor(negmask,:) = repmat([1 0 0], sum(negmask), 1);
%% edges incident on end nodes
fedges = find(T(:,fend) + F(:,fend));
fedges = fedges(fedges ~= bid);
tedges = find(T(:,tend) + F(:,tend));
tedges = tedges(tedges ~= bid);

%% neighbors of end nodes
fneighbors = find(Adj(fend,:)).';
fneighbors = fneighbors(fneighbors ~= tend);
tneighbors = find(Adj(tend,:)).';
tneighbors = tneighbors(tneighbors ~= fend);

%% compare branch parameters
r =  Sb.g./(Sb.g.^2 + Sb.b.^2);
x = -Sb.b./(Sb.g.^2 + Sb.b.^2);
bcharge = F'*Sb.bc/2 + T'*Sb.bc/2;

nodeids = [fend; tend; fneighbors; tneighbors];
edgeids = [bid; fedges; tedges];

pqmask  = bidx.pq(nodeids);
pvmask  = bidx.pv(nodeids);
rfmask  = bidx.ref(nodeids);
%% graph of neighborhood
tmpmap  = sparse(nodeids, 1, 1:length(nodeids));
% tmpadj  = sparse(full(tmpmap([flist(edgeids);tlist(edgeids)])), ...
%                  full(tmpmap([tlist(edgeids);flist(edgeids)])), 1);
tmpadj  = sparse(full(tmpmap(flist(edgeids))), full(tmpmap(tlist(edgeids))), 1, length(nodeids), length(nodeids));
Grph = digraph(tmpadj>0,cellstr(num2str(nodeids)));
edgelabels = cell(size(Grph.Edges,1),1);
for k = 1:length(edgeids)
    tmpid = Grph.findedge(full(tmpmap(flist(edgeids(k)))), full(tmpmap(tlist(edgeids(k)))));
    if isempty(edgelabels{tmpid})
        edgelabels{tmpid} = num2str(edgeids(k));
    else
        edgelabels{tmpid} = sprintf('%s,%d',edgelabels{tmpid},edgeids(k));
    end
end
figure;
h = plot(Grph,'edgelabel',edgelabels);
highlight(h,full(tmpmap(fend)),full(tmpmap(tend)))
highlight(h,full(tmpmap(fend)),full(tmpmap(tend)), 'EdgeColor', 'r')
highlight(h,full(tmpmap(nodeids(pqmask))))
highlight(h,full(tmpmap(nodeids(pvmask))), 'NodeColor', 'g', 'MarkerSize',10)
highlight(h,full(tmpmap(nodeids(rfmask))), 'NodeColor', 'k')
%%
figure;
cmap = colormap('lines');
%%%% resistance
subplot(5,2,1)
bar(1,r(bid),'FaceColor', cmap(1,:));
hold on;
bar(2:2+length(fedges)-1,r(fedges), 'FaceColor', cmap(2,:))
bar(2+length(fedges):2+length(fedges)+ length(tedges)-1,r(tedges),'FaceColor', cmap(3,:))
ax = gca;
ax.XTick = 1:length(edgeids);
ax.XTickLabel = cellstr(num2str(edgeids)).';
ylabel('resistance p.u')
legend('bid', 'fedges', 'tedges')

%%%% reactance
subplot(5,2,3)
bar(1,x(bid),'FaceColor', cmap(1,:));
hold on;
bar(2:2+length(fedges)-1,x(fedges),'FaceColor', cmap(2,:))
bar(2+length(fedges):2+length(fedges)+ length(tedges)-1,x(tedges),'FaceColor', cmap(3,:))
ax = gca;
ax.XTick = 1:length(edgeids);
ax.XTickLabel = cellstr(num2str(edgeids)).';
ylabel('reactance p.u')


%%%% r/x
subplot(5,2,5)
bar(1,r(bid)./x(bid),'FaceColor', cmap(1,:));
hold on;
bar(2:2+length(fedges)-1,r(fedges)./x(fedges),'FaceColor', cmap(2,:))
bar(2+length(fedges):2+length(fedges)+ length(tedges)-1,r(tedges)./x(tedges),'FaceColor', cmap(3,:))
ax = gca;
ax.XTick = 1:length(edgeids);
ax.XTickLabel = cellstr(num2str(edgeids)).';
ylabel('r/x')


%%%% tau
subplot(5,2,7)
bar(1,Sb.tau(bid),'FaceColor', cmap(1,:));
hold on;
bar(2:2+length(fedges)-1,Sb.tau(fedges), 'FaceColor', cmap(2,:))
bar(2+length(fedges):2+length(fedges)+ length(tedges)-1,Sb.tau(tedges),'FaceColor', cmap(3,:))
ax = gca;
ax.XTick = 1:length(edgeids);
ax.XTickLabel = cellstr(num2str(edgeids)).';
ylabel('\tau')

%%%% tshift
subplot(5,2,9)
bar(1,Sb.tshift(bid),'FaceColor', cmap(1,:));
hold on;
bar(2:2+length(fedges)-1,Sb.tshift(fedges), 'FaceColor', cmap(2,:))
bar(2+length(fedges):2+length(fedges)+ length(tedges)-1,Sb.tshift(tedges),'FaceColor', cmap(3,:))
ax = gca;
ax.XTick = 1:length(edgeids);
ax.XTickLabel = cellstr(num2str(edgeids)).';
xlabel('edge id')
ylabel('\delta')


% node properties

%%%% shunt reactiance (line charging)
subplot(5,2,2)
bar(1, bcharge(fend), 'FaceColor', cmap(1,:))
hold on;
bar(2, bcharge(tend), 'FaceColor', cmap(2,:))
bar(3:3+length(fneighbors)-1, bcharge(fneighbors), 'FaceColor', cmap(3,:))
bar(3+length(fneighbors):3+length(fneighbors)+length(tneighbors)-1,...
    bcharge(tneighbors), 'FaceColor', cmap(4,:))
ax = gca;
ax.XTick = 1:length(nodeids);
ax.XTickLabel = cellstr(num2str(nodeids));
ylabel('sum(bc/2)')
legend('fend', 'tend', 'fneighbors', 'tneighbors')

%%%% shunt reactiance (not line charging)
subplot(5,2,4)
bar(1, Sb.bsh(fend), 'FaceColor', cmap(1,:))
hold on;
bar(2, Sb.bsh(tend), 'FaceColor', cmap(2,:))
bar(3:3+length(fneighbors)-1, Sb.bsh(fneighbors), 'FaceColor', cmap(3,:))
bar(3+length(fneighbors):3+length(fneighbors)+length(tneighbors)-1,...
    Sb.bsh(tneighbors), 'FaceColor', cmap(4,:))
ax = gca;
ax.XTick = 1:length(nodeids);
ax.XTickLabel = cellstr(num2str(nodeids));
ylabel('Bsh')
legend('fend', 'tend', 'fneighbors', 'tneighbors')

%%%% shunt conductance
subplot(5,2,6)
bar(1, Sb.gsh(fend), 'FaceColor', cmap(1,:))
hold on;
bar(2, Sb.gsh(tend), 'FaceColor', cmap(2,:))
bar(3:3+length(fneighbors)-1, Sb.gsh(fneighbors), 'FaceColor', cmap(3,:))
bar(3+length(fneighbors):3+length(fneighbors)+length(tneighbors)-1,...
    Sb.gsh(tneighbors), 'FaceColor', cmap(4,:))
ax = gca;
ax.XTick = 1:length(nodeids);
ax.XTickLabel = cellstr(num2str(nodeids));
ylabel('Gsh')

%%%% real load
subplot(5,2,8)
bar(1, Sp.Pd(fend), 'FaceColor', cmap(1,:))
hold on;
bar(2, Sp.Pd(tend), 'FaceColor', cmap(2,:))
bar(3:3+length(fneighbors)-1, Sp.Pd(fneighbors), 'FaceColor', cmap(3,:))
bar(3+length(fneighbors):3+length(fneighbors)+length(tneighbors)-1,...
    Sp.Pd(tneighbors), 'FaceColor', cmap(4,:))
ax = gca;
ax.XTick = 1:length(nodeids);
ax.XTickLabel = cellstr(num2str(nodeids));
ylabel('Pd')
xlabel('node id')

%%%% reactive load
subplot(5,2,10)
bar(1, Sp.Qd(fend), 'FaceColor', cmap(1,:))
hold on;
bar(2, Sp.Qd(tend), 'FaceColor', cmap(2,:))
bar(3:3+length(fneighbors)-1, Sp.Qd(fneighbors), 'FaceColor', cmap(3,:))
bar(3+length(fneighbors):3+length(fneighbors)+length(tneighbors)-1,...
    Sp.Qd(tneighbors), 'FaceColor', cmap(4,:))
ax = gca;
ax.XTick = 1:length(nodeids);
ax.XTickLabel = cellstr(num2str(nodeids));
ylabel('Qd')
xlabel('node id')