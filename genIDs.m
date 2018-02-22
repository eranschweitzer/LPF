function IDs = genIDs(varargin)
% generate a matrix of all the iteration IDS
iters = varargin_parse(varargin,'iterID',0:3);

IDs = int8(zeros(4*(5*2)*(8*2)*(4*2)*(3*2)*(5*2)*(3*2), 7));
nv = [2,4,5,7];
i = 1;
for iterID = iters
    for Art= [0,1,3] %0:4
        for Artb = 0:1
            for Aru = [0,1,3,5,6]%0:7
                for Arub = 0:1
                    for Amt = [0,1,5]%[0,1,2,5]
                        for Amtb = 0:1
                            for Amu = [0,1]%0:2
                                for Amub = 0:1
                                    for br = [0,1,3]%0:4
                                        for brb = 0:1
                                            for bm = [0,1]%0:2
                                                for bmb = 0:1 
                                                    % cases to skip
                                                    if (Aru > 4) && (Arub == 1)
                                                        % cases above 4 do not have b in them
                                                        continue
                                                    elseif (Amt > 0) && (Amtb == 1)
                                                        % cases other than 0 do not have b in them
                                                        continue
                                                    elseif (iterID > 1) && ismember(Art,nv) && ismember(Aru,nv) && ismember(Amt,nv) && ismember(Amu,nv) && ismember(br,nv) && ismember(bm,nv)
                                                        % no u0 iteration when u0 is neglected in all parts
                                                        continue
                                                    else
                                                        IDs(i,:) = [iterID, Art*10+Artb, Aru*10+Arub, Amt*10+Amtb, Amu*10+Amub, br*10+brb, bm*10+bmb];
                                                        i = i + 1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
IDs = IDs(1:i-1,:);