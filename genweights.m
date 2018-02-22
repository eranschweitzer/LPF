function Gw = genweights(x,mode)


if strcmp(mode,'default')
    Gw = 1;
elseif isscalar(mode)
    Gw = diags(1./sqrt((1-(x./mode).^2).^2 + (2*0.7*(x./mode)).^2));
elseif contains(mode,'avg')
    if contains(mode,'avgstd')
        s = str2double(mode(7:end));
    else
        s = 0;
    end
    xn = mean(x) + s*std(x);
    Gw = diags(1./sqrt((1-(x./xn).^2).^2 + (2*0.7*(x./xn)).^2));
else
    error('unknown weight method')
end
