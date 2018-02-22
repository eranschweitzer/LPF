function T = resprep_flow(ids, flids, Cf, caseid)

f = fieldnames(Cf.pf{1});
N = 2*length(f)*(length(Cf.pf) + length(Cf.qf));
idcol = cell(N,1);
prop  = cell(N,1);
crit  = cell(N,1);
val   = zeros(N,1);

idx = 1;
for a = {'p','q'}
	for l = 1:length(flids.(a{1}))
		for b = {'f', 't'}
			for k = 1:length(f)
				prop{idx} = [a{1}, b{1}];
				crit{idx} = f{k};
				val(idx)  = Cf.(prop{idx}){l}.(f{k});
				idcol{idx}= flids.(a{1}){l};
				idx = idx + 1;
			end
		end
	end
end

casecol = caseid*ones(N,1);
pfidcol = char(ones(N,1)*ids2str(ids));

T = table(casecol, pfidcol, idcol, prop, crit, val, 'VariableNames', {'caseid', 'pfid', 'id', 'prop', 'criteria', 'value'});

