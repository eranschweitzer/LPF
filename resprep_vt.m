function T = resprep_vt(ids,Cv, Ct, caseid)

f = fieldnames(Cv);
lf = length(f);
prop = cell(lf*2,1);
crit = cell(lf*2,1);
val  = zeros(lf*2,1);

for k = 1:lf
	prop{k} = 'vol';
	crit{k} = f{k};
	val(k)  = Cv.(f{k});
end

for k = 1:lf
	prop{k + lf} = 'ang';
	crit{k + lf} = f{k};
	val(k  + lf) = Ct.(f{k});
end

casecol = caseid*ones(lf*2,1);
idcol    = char(ones(lf*2,1)*ids2str(ids));

T = table(casecol, idcol, prop, crit, val, 'VariableNames', {'caseid', 'id', 'prop', 'criteria', 'value'});
