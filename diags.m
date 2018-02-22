function X = diags(x)
% turns vector x into a sparse diagonal matrix X

X = sparse(1:length(x),1:length(x),x);