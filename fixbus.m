function [bidx,N,matI,matU] = fixbus(busmask,u0,bidx,N)


%% fix u0 voltages and change pv and pq indecies
bidx.pq(busmask) = false;
bidx.pv(busmask) = true;
N.pv = sum(bidx.pv); %update count of PV buses
N.pq = sum(bidx.pq);
%% matrices fixing PV and Ref quantities
[matI,matU] = pv_ref_mats(bidx,N,u0);