function convg = linpf2ac_convergence(mpc,vars,bidx)
%%% check whether the ac powerflow converges when the values from the
%%% linear powerflows are used. Use both the dc and the linpf.

%% basic setup
define_constants;
mpopt   = mpoption; mpopt.out.all = 0;
%% DC powerflow
if nargin > 2
    mpcdc = mpc;
    % flat start conditions
    mpcdc.bus(bidx.pq,VM) = 1;
    mpcdc.bus(bidx.pq | bidx.pv,VA) = 0;
    mpcdc = rundcpf(mpcdc,mpopt); % 1) solve DC
    mpcdc = runpf(mpcdc,mpopt);   % 2) solve ACPF
end
%% linpf
mpctest = mpc;
% linpf solution
mpctest.bus(:,VM) = vars.v;
mpctest.bus(:,VA) = vars.theta*180/pi;
mpctest = runpf(mpctest,mpopt);

%% results
if nargin > 2
    convg = struct('dc', mpcdc.success, 'linpf', mpctest.success);
else
    convg = struct('linpf',mpctest.success);
end