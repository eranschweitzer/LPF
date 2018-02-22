function T = result2table(colnames,casename,ids,convg,Cv,Ct,Cpf,Cpt,Cqf,Cqt,pfids)

if nargin == 10
    pfids = 0:5;
end   
%% setup structures
c = {'cor','max','avg','del'};
% N = 0+4+4+4*numel(Cpf)*4;
N = 4+4+4*length(pfids)*2*4;
casecol  = char(ones(N,1)*casename);
idcol    = char(ones(N,1)*ids2str(ids));
prop     = cell(N,1);
criteria = cell(N,1);
val      = zeros(N,1);
%% convergence
% prop{1}     = num2str(ids.iter);
% criteria{1} = 'con';
% val(1)      = convg;
%% voltage data
ptr = 1;
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
for ID = (pfids+1)%1:size(Cpf,1)
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
for ID = (pfids+1)%1:size(Cpt,1)
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
for ID = (pfids+1)%1:size(Cqf,1)
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
for ID = (pfids+1)%1:size(Cqt,1)
    for btype = 1:2
        for i = 1:length(c)
            prop{ptr}    = ['qt', num2str(ID-1), num2str(btype-1)];
            criteria{ptr}= propmap(c{i});
            val(ptr)     = Cqt{ID,btype}.(c{i});
            ptr = ptr + 1;
        end
    end
end

%% create table
T = table(casecol,idcol,prop,criteria,val,'VariableNames',colnames);
