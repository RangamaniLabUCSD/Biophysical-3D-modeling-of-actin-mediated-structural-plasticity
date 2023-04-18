function dvar_dt = f_var(t,var,fD,dd,fU,J,P)

%     [F_d,A,C]
    fd = reshape(var(1:P.n_x3),P.n_x,P.n_x,P.n_x);
    a = reshape(var((1*P.n_x3+1):2*P.n_x3),P.n_x,P.n_x,P.n_x);
    c = reshape(var((2*P.n_x3+1):3*P.n_x3),P.n_x,P.n_x,P.n_x);
    
    fD_x = reshape(fD(:,1),P.n_x,P.n_x,P.n_x);
    fD_y = reshape(fD(:,2),P.n_x,P.n_x,P.n_x);
    fD_z = reshape(fD(:,3),P.n_x,P.n_x,P.n_x);
    
    aux_fD_z = zeros(P.n_x,P.n_x,P.n_x);
    aux_fD_z(:,:,1) = fD_z(:,:,2)-fD_z(:,:,end);
    aux_fD_z(:,:,2:end-1) = fD_z(:,:,3:end)- fD_z(:,:,1:end-2);
    aux_fD_z(:,:,end) = fD_z(:,:,1)-fD_z(:,:,end-1);
    
    fU_x = reshape(fU(:,1),P.n_x,P.n_x,P.n_x);
    fU_y = reshape(fU(:,2),P.n_x,P.n_x,P.n_x);
    fU_z = reshape(fU(:,3),P.n_x,P.n_x,P.n_x);
    
    aux_fU_z = zeros(P.n_x,P.n_x,P.n_x);
    aux_fU_z(:,:,1) = fU_z(:,:,2)-fU_z(:,:,end);
    aux_fU_z(:,:,2:end-1) = fU_z(:,:,3:end)- fU_z(:,:,1:end-2);
    aux_fU_z(:,:,end) = fU_z(:,:,1)-fU_z(:,:,end-1);
    
    %Basal influx
    aux_Omega = (dd(:,1)<-P.d_th & P.Z(:) > P.head_ini);
    aux_Omega = aux_Omega.*P.B;
    Omega = reshape(aux_Omega,P.n_x,P.n_x,P.n_x).*sum(P.B)./sum(aux_Omega);
%     stimulus-triggered influx for barbed ends
    aux_Omega2 = ((dd(:,1)<-P.d_th & dd(:,1)>-P.d_th5) & (P.Z(:)>P.head_stim & P.Z(:)<P.head_stim2));
    Omega2 = reshape(aux_Omega2,P.n_x,P.n_x,P.n_x).*sum(P.B2)./sum(aux_Omega2);
%     stimulus-triggered influx for Arp2/3
    aux_Omega3 = ((dd(:,1)<-P.d_th & dd(:,1)>-P.d_th2) & (P.Z(:)>P.head_stim & P.Z(:)<P.head_stim2));
    Omega3 = reshape(aux_Omega3,P.n_x,P.n_x,P.n_x).*sum(P.B3)./sum(aux_Omega3);
%     stimulus-triggered influx for cofilin  
    aux_Omega4 = ((dd(:,1)<-P.d_th3 & dd(:,1)>-P.d_th4) & (P.Z(:)>P.head_stim & P.Z(:)<P.head_stim2));
    Omega4 = reshape(aux_Omega4,P.n_x,P.n_x,P.n_x).*sum(P.B4)./sum(aux_Omega4);

    f_nuc = P.k_nuc*P.psi_1*a.*fd;
    f_sev = P.k_sev*P.psi_1*fd.*(c.^P.n)./(P.k_n + c.^P.n);
% calculate rho
    aux_vel = -(P.nu/(2*P.delta_x)).*...
        (fd.*[fD_x(:,2,:)-fD_x(:,end,:),...
        fD_x(:,3:end,:)-fD_x(:,1:end-2,:),...
        fD_x(:,1,:)-fD_x(:,end-1,:)] + ...
        fd.*[fD_y(2,:,:)-fD_y(end,:,:);
        fD_y(3:end,:,:)-fD_y(1:end-2,:,:);
        fD_y(1,:,:)-fD_y(end-1,:,:)]+...
        fd.*aux_fD_z)...
        + (P.a_rp/(2*P.delta_x)).*...
        (fd.*[fU_x(:,2,:)-fU_x(:,end,:),...
        fU_x(:,3:end,:)-fU_x(:,1:end-2,:),...
        fU_x(:,1,:)-fU_x(:,end-1,:)] + ...
        fd.*[fU_y(2,:,:)-fU_y(end,:,:);
        fU_y(3:end,:,:)-fU_y(1:end-2,:,:);
        fU_y(1,:,:)-fU_y(end-1,:,:)]+...
        fd.*aux_fU_z);
    
    
    dfd_dt = aux_vel + P.psi_0*(f_nuc + f_sev + P.I_B*Omega+ J(1)*Omega2) - P.k_B*fd;
    
    da_dt = P.zeta_b*P.psi_1*a.*aux_vel - P.k_A*a + P.I_A*Omega - f_nuc +J(2)*Omega3;
    
    dc_dt =  P.zeta_b*P.psi_1*c.*aux_vel - P.k_C*c + P.I_C*Omega - f_sev + J(3)*Omega4; 
        
    dvar_dt = [dfd_dt(:);da_dt(:);dc_dt(:)];
end