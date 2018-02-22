function ids = flowids()
p = cell(18,1); 
q = cell(34,1);
pidx = 1;
qidx = 1;
for opt = 1:2
	for AB = 1:2
		for simp = 0:1
			for lin = 0:1
				for g = 0:1
					q{qidx} = sprintf('%d%d%d%d%d', opt, AB, simp, lin, g);	
					qidx = qidx + 1;
					if opt == 2
						continue
					end
					p{pidx} = sprintf('%d%d%d%d%d', opt, AB, simp, lin, g);	
					pidx = pidx + 1;
				end
			end
		end
	end
end

p{pidx} = '30000';
q{qidx} = '30000';
pidx = pidx + 1;
qidx = qidx + 1;
p{pidx} = '40000';
q{qidx} = '40000';

ids = struct('p', {p}, 'q', {q});
