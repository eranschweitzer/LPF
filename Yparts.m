function [Yff,Yft,Ytf,Ytt,Ysh] = Yparts(mpc)
% return the parts of the admittance matrix. Essentially copied from
% Matpower's makeYbus

define_constants;
nl = size(mpc.branch,1);
% for each branch, compute the elements of the branch admittance matrix where
%
%      | If |   | Yff  Yft |   | Vf |
%      |    | = |          | * |    |
%      | It |   | Ytf  Ytt |   | Vt |
%
stat = mpc.branch(:, BR_STATUS);                    %% ones at in-service branches
Ys = stat ./ (mpc.branch(:, BR_R) + 1i * mpc.branch(:, BR_X));  %% series admittance
Bc = stat .* mpc.branch(:, BR_B);                           %% line charging susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
t = find(mpc.branch(:, TAP));                       %% indices of non-zero tap ratios
tap(t) = mpc.branch(t, TAP);                        %% assign non-zero tap ratios
tap = tap .* exp(1i*pi/180 * mpc.branch(:, SHIFT)); %% add phase shifters
Ytt = Ys + 1i*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

% compute shunt admittance
% if Psh is the real power consumed by the shunt at V = 1.0 p.u.
% and Qsh is the reactive power injected by the shunt at V = 1.0 p.u.
% then Psh - j Qsh = V * conj(Ysh * V) = conj(Ysh) = Gs - j Bs,
% i.e. Ysh = Psh + j Qsh, so ...
Ysh = (mpc.bus(:, GS) + 1i * mpc.bus(:, BS)) / mpc.baseMVA; %% vector of shunt admittances