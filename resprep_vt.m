function T = resprep_vt(ids,Cv, Ct, Cdt, caseid)

f = fieldnames(Cv);
lf = length(f);
prop = cell(lf*3,1);
crit = cell(lf*3,1);
val  = zeros(lf*3,1);

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

for k = 1:lf
	prop{k + 2*lf} = 'del';
	crit{k + 2*lf} = f{k};
	val(k  + 2*lf) = Cdt.(f{k});
end

casecol = caseid*ones(lf*3,1);
idcol    = char(ones(lf*3,1)*ids2str(ids));

T = table(casecol, idcol, prop, crit, val, 'VariableNames', {'caseid', 'id', 'prop', 'criteria', 'value'});
