function [num,den] = diagpade(p)
%%% p are the polynomial coefficients in decreasing order

deg = (length(p) - 1)/2;

coeffmat = cell(length(p)-1 - deg,1);
for k = 1:length(p)-(deg+1)
    coeffmat{k} = p(end-k:-1:end-(deg+k-1));
end

den = [cell2mat(coeffmat)\-p(end-(deg+1):-1:1).'; 1];

coeffmat = cell(1,deg+1);
for k = 0:length(coeffmat)-1
    coeffmat{k+1} = diag(ones(deg+1 - k,1),-k)*p(end:-1:end-deg)';
end

num = flipud(cell2mat(coeffmat)*flipud(den));