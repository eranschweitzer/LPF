function result2db_old(tablename,casename,ids,convg,Cv,Ct,Cpf,Cpt,Cqf,Cqt)

%% create connection
% since we are adding one row at a time, we leave the autocommit setting on
conn = dbconn();

%% setup structures
c = {'corr','max','avg','delta'};
N = 1+4+4+4*numel(Cpf)*4;
% colnames = cell(1,3+4+4+4*numel(Cpf)*4);
% data     = cell(1,3+4+4+4*numel(Cpf)*4);
% colnames{1} = 'casename';  data{1} = casename;
% colnames{2} = 'id';    data{2} = ids2str(ids);
% colnames{3} = 'convg'; data{3} = logical(convg);
colnames = {'casename', 'id', 'prop', 'criteria', 'value'};
casecol  = char(ones(N,1)*casename);
idcol    = char(ones(N,1)*ids2str(ids));
prop     = cell(N,1);
criteria = cell(N,1);
val      = zeros(N,1);
%% convergence
prop{1}     = num2str(ids.iter);
criteria{1} = 'con';
val(1)      = convg;
ptr = 2;
%% voltage data
for i = 1:length(c)
    prop{ptr}    = 'volt';
    criteria{ptr}= propmap(c{i});
    val(ptr)     = Cv.(c{i});
    ptr = ptr + 1;
end

%% angle data
for i = 1:length(c)
    prop{ptr}    = 'ang';
    criteria{ptr}= propmap(c{i});
    val(ptr)     = Ct.(c{i});
    ptr = ptr + 1;
end

%% Pf data
for ID = 1:size(Cpf,1)
    for btype = 1:2
        for i = 1:length(c)
            prop{ptr}    = ['pf', num2str(ID-1), num2str(btype-1)];
            criteria{ptr}= propmap(c{i});
            val(ptr)     = Cpf{ID,btype}.(c{i});
            ptr = ptr + 1;
        end
    end
end

%% Pt data
for ID = 1:size(Cpt,1)
    for btype = 1:2
        for i = 1:length(c)
            prop{ptr}    = ['pt', num2str(ID-1), num2str(btype-1)];
            criteria{ptr}= propmap(c{i});
            val(ptr)     = Cpt{ID,btype}.(c{i});
            ptr = ptr + 1;
        end
    end
end

%% Qf data
for ID = 1:size(Cqf,1)
    for btype = 1:2
        for i = 1:length(c)
            prop{ptr}    = ['qf', num2str(ID-1), num2str(btype-1)];
            criteria{ptr}= propmap(c{i});
            val(ptr)     = Cqf{ID,btype}.(c{i});
            ptr = ptr + 1;
        end
    end
end

%% Qt data
for ID = 1:size(Cqt,1)
    for btype = 1:2
        for i = 1:length(c)
            prop{ptr}    = ['qt', num2str(ID-1), num2str(btype-1)];
            criteria{ptr}= propmap(c{i});
            val(ptr)     = Cqt{ID,btype}.(c{i});
            ptr = ptr + 1;
        end
    end
end

%% insert data
datainsert(conn,tablename,colnames,table(casecol,idcol,prop,criteria,val,'VariableNames',colnames));

%% close connection
close(conn)