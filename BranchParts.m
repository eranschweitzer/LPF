function S = BranchParts(mpc)
% return the parts of the admittance matrix. Essentially copied from
% Matpower's makeYbus
S = struct();
define_constants;
nl = size(mpc.branch,1);
% for each branch, compute the elements of the branch admittance matrix where
%
%      | If |   | Yff  Yft |   | Vf |
%      |    | = |          | * |    |
%      | It |   | Ytf  Ytt |   | Vt |
%
stat = mpc.branch(:, BR_STATUS);                    %% ones at in-service branches
S.g  = ( stat.*mpc.branch(:,BR_R))./(mpc.branch(:, BR_R).^2 + mpc.branch(:, BR_X).^2);
S.b  = (-stat.*mpc.branch(:,BR_X))./(mpc.branch(:, BR_R).^2 + mpc.branch(:, BR_X).^2);
S.balt = -1./mpc.branch(:,BR_X);
S.bc = stat .* mpc.branch(:, BR_B);                           %% line charging susceptance
S.tau = ones(nl, 1);                              %% default tap ratio = 1
t = find(mpc.branch(:, TAP));                       %% indices of non-zero tap ratios
S.tau(t) = mpc.branch(t, TAP);                        %% assign non-zero tap ratios
S.tshift = pi/180 * mpc.branch(:, SHIFT); %% phase shifters

% compute shunt admittance
% if Psh is the real power consumed by the shunt at V = 1.0 p.u.
% and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
% then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
% i.e. Ysh = Psh + j Qsh, so ...
% Ysh = (mpc.bus(:, GS) + 1i * mpc.bus(:, BS)) / mpc.baseMVA; %% vector of shunt admittances
S.gsh = mpc.bus(:, GS) / mpc.baseMVA; %% vector of shunt conductance
S.bsh = mpc.bus(:, BS) / mpc.baseMVA; %% vector of shunt susceptance