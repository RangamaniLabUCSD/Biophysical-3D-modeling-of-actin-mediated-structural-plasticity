function [S_int,var_0] = solve_system_3D(S,T_s,var_0,t,delta_t,J,v_mm,P)
    dd = dist_r_s_3D(S,T_s,P);
    [fu,fd] = f_u(S,T_s,dd,P);
    t_span = [t t+delta_t];
    [~,vars] = ode45(@(t,var) f_var(t,var,fd,dd,fu,J,v_mm,P), t_span, var_0);
    var_0 = vars(end,:);
   
    f_mem = force_membrane_3D(S,T_s,P);
    f_actin = force_actin_3D(S,T_s,dd,var_0(1:P.n_x3)',P.alpha,P.beta,P);
    S_int = S + delta_t*P.zeta*(f_mem+f_actin);%P.zeta = 1/omega
  
end