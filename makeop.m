function op = makeop(ids, u0)

if ~iscell(ids)
    ids = {ids};
end
if ~iscell(u0)
    u0 = {u0};
end

op = cell(numel(ids)*numel(u0),1);
ptr = 1
for ki = 1:length(ids)
    for ku = 1:length(u0)
        op{ptr} = struct('ids', ids{ki}, 'u0', u0{ku});
        ptr = ptr + 1;
    end
end

