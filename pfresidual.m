function err = pfresidual(v, Ybus, Pinj, Qinj)
%%% err = pfresidual(v, Ybus, Pinj, Qinj;

err = diags(conj(v))*Ybus*v - (Pinj - 1i*Qinj);