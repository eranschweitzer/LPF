function S = eval_criteria(x,xtrue)

S = struct();
S.cor      = corr(x,xtrue);
[S.max,idx] = max(abs(x-xtrue));
S.avg       = mean(abs(x-xtrue));
S.del     = 100*abs((x(idx) - xtrue(idx))/xtrue(idx));
S.rms       = 1/sqrt(length(x)) * norm(x-xtrue);
