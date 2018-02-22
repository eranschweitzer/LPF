function uavg = ufix(u0,E,busmask)
% sum the neighbor u0 of violating u0. 
% exclude the violating u0 from the sum.
% divide by degree (excluding violating buses) to get the average

A = -(E'*E);
A = A - diags(diag(A));
A(A>0) = 1;

uavg = (A(busmask,~busmask)*u0(~busmask))./sum(A(busmask,~busmask),2);
uavg(isnan(uavg) | isinf(uavg)) = 0;%mean(u0(~busmask));