function expansion_error(vtrue,u0)

vtrue.v - exp(u0).*(1 + log(vtrue.v) - u0)
(vtrue.v).^(-1) - exp(-u0).*(1 - (log(vtrue.v) - u0))

(vtrue.v).^(-2) - exp(-2*u0).*(1 - 2*(log(vtrue.v) - u0))

T*vtrue.v.*exp(-(F-T)*vtrue.t) - exp(T*u0).*(1 + T*(log(vtrue.v) - u0) - 1i*(F-T)*vtrue.t)

T*vtrue.v.*exp(-(F-T)*vtrue.t) - exp(T*u0).*(1 + T*(log(vtrue.v) - u0) - 1i*(F-T)*vtrue.t - 0.5*((F-T)*vtrue.t).^2)

a = exp(1i*(F-T)*vtrue.t);
b = (1 + 1i*(F-T)*vtrue.t - 0.5*((F-T)*vtrue.t).^2 -1i*((F-T)*vtrue.t).^3/6 + ((F-T)*vtrue.t).^4/24);
