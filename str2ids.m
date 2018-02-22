function ids = str2ids(s)

if length(s) == 13
    ids.iter = str2double(s(1));
    ids.Art  = str2double(s(2));
    ids.Artb = str2double(s(3));
    ids.Aru  = str2double(s(4));
    ids.Arub = str2double(s(5));
    ids.Amt  = str2double(s(6));
    ids.Amtb = str2double(s(7));
    ids.Amu  = str2double(s(8));
    ids.Amub = str2double(s(9));
    ids.br   = str2double(s(10));
    ids.brb  = str2double(s(11));
    ids.bm   = str2double(s(12));
    ids.bmb  = str2double(s(13));
elseif length(s) == 12
    fields = {'iter', 'Q', 'AB', 'c', 'gr', 'gi',...
              'Omr', 'Gami', 'Kr', 'Omi', 'Gamr', 'Ki'};
    ids = struct();
    for k = 1:length(fields)
        ids.(fields{k}) = str2double(s(k));
    end
elseif length(s) == 5
    fields = {'Q', 'AB', 'simp', 'lin', 'g'};
    ids = struct();
    for k = 1:length(fields)
        ids.(fields{k}) = str2double(s(k));
    end
else
    error('Incorrect ID string. Number of elements is %d',  length(s))
end
