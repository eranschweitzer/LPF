function S = unpackIDs(ID)

if ischar(ID)
	S = str2ids(ID);
elseif iscell(ID)
	if numel(ID) == 1
		S = str2ids(ID{1});
	else
		error('Too many cell elements %d', numel(ID))
	end
else
	S = struct();
	S.iter = ID(1);
	[S.Art, S.Artb] = btype_extract(ID(2));
	[S.Aru, S.Arub] = btype_extract(ID(3));
	[S.Amt, S.Amtb] = btype_extract(ID(4));
	[S.Amu, S.Amub] = btype_extract(ID(5));
	[S.br,  S.brb]  = btype_extract(ID(6));
	[S.bm,  S.bmb]  = btype_extract(ID(7));
end
