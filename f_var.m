function dvar_dt = f_var(t,var,fD,dd,fU,J,v_mm,P)

%     [B,a,c] B=barbed ends (n^3,1), a=Arp2,3
    
    B = reshape(var(1:P.n_x3),P.n_x,P.n_x,P.n_x);
    a = reshape(var((1*P.n_x3+1):2*P.n_x3),P.n_x,P.n_x,P.n_x);
    c = reshape(var((2*P.n_x3+1):3*P.n_x3),P.n_x,P.n_x,P.n_x);
    
%     absorbing boudnary conditions
    aux_poly =  (P.Z(:)>P.head_stim+P.delta_x);
       
    aux_Bound = (dd(:,1)<-P.delta_x);
    aux_Bound = reshape(aux_poly,P.n_x,P.n_x,P.n_x).*reshape(aux_Bound,P.n_x,P.n_x,P.n_x);%
    
    fD_x = reshape(fD(:,1),P.n_x,P.n_x,P.n_x);
    fD_y = reshape(fD(:,2),P.n_x,P.n_x,P.n_x);
    fD_z = reshape(fD(:,3),P.n_x,P.n_x,P.n_x);
    
    fD_x_B = B.*fD_x;
    fD_y_B = B.*fD_y;
    fD_z_B = B.*fD_z;
    
    fD_x_a = a.*fD_x;
    fD_y_a = a.*fD_y;
    fD_z_a = a.*fD_z;
    
    fD_x_c = c.*fD_x;
    fD_y_c = c.*fD_y;
    fD_z_c = c.*fD_z;

    
    aux_fD_z_B = zeros(P.n_x,P.n_x,P.n_x);

    aux_fD_z_B(:,:,1) = fD_z_B(:,:,2).*aux_Bound(:,:,1);
    aux_fD_z_B(:,:,2:end-1) = (fD_z_B(:,:,3:end)- fD_z_B(:,:,1:end-2)).*aux_Bound(:,:,2:end-1);
    aux_fD_z_B(:,:,end) = -fD_z_B(:,:,end-1).*aux_Bound(:,:,end);
    
    aux_fD_z_c = zeros(P.n_x,P.n_x,P.n_x);
    aux_fD_z_c(:,:,1) = fD_z_c(:,:,2).*aux_Bound(:,:,1);
    aux_fD_z_c(:,:,2:end-1) = (fD_z_c(:,:,3:end)- fD_z_c(:,:,1:end-2)).*aux_Bound(:,:,2:end-1);
    aux_fD_z_c(:,:,end) = -fD_z_c(:,:,end-1).*aux_Bound(:,:,end);
    
    aux_fD_z_a = zeros(P.n_x,P.n_x,P.n_x);
    aux_fD_z_a(:,:,1) = fD_z_a(:,:,2).*aux_Bound(:,:,1);
    aux_fD_z_a(:,:,2:end-1) = (fD_z_a(:,:,3:end)- fD_z_a(:,:,1:end-2)).*aux_Bound(:,:,2:end-1);
    aux_fD_z_a(:,:,end) = -fD_z_a(:,:,end-1).*aux_Bound(:,:,end);
    
    fU_x = reshape(fU(:,1),P.n_x,P.n_x,P.n_x);
    fU_y = reshape(fU(:,2),P.n_x,P.n_x,P.n_x);
    fU_z = reshape(fU(:,3),P.n_x,P.n_x,P.n_x);
    
    fU_x_B = B.*fU_x;
    fU_y_B = B.*fU_y;
    fU_z_B = B.*fU_z;
    

    aux_fU_z_B = zeros(P.n_x,P.n_x,P.n_x);
    aux_fU_z_B(:,:,1) = fU_z_B(:,:,2).*aux_Bound(:,:,1);
    aux_fU_z_B(:,:,2:end-1) = (fU_z_B(:,:,3:end)- fU_z_B(:,:,1:end-2)).*aux_Bound(:,:,2:end-1);
    aux_fU_z_B(:,:,end) = -fU_z_B(:,:,end-1).*aux_Bound(:,:,end);
    
   

    %Basal influx
  
    %Basal influx
    Omega = reshape(P.B,P.n_x,P.n_x,P.n_x);
%     stimulus-triggered influx for barbed ends
    aux_Omega2 = ((dd(:,1)<-P.d_th) & (P.Z(:)>P.head_stim & P.Z(:)<P.head_stim2));
    Omega2 = reshape(aux_Omega2,P.n_x,P.n_x,P.n_x).*sum(P.B2)./sum(aux_Omega2);
    

    f_nuc = P.k_nuc*P.psi_1*a.*B;
    f_sev = P.k_sev*P.psi_1*B.*(c.^P.n)./(P.k_n + c.^P.n);
% calculate rho
    nu = P.nu*(J(1)+P.I_B)/P.I_B;
    aux_vel_B = -(nu/(2*P.delta_x)).*...
            ([fD_x_B(:,2,:).*aux_Bound(:,1,:),...
            (fD_x_B(:,3:end,:)-fD_x_B(:,1:end-2,:)).*aux_Bound(:,2:end-1,:),...
            -fD_x_B(:,end-1,:).*aux_Bound(:,end,:)] + ...
            [fD_y_B(2,:,:).*aux_Bound(1,:,:);
            (fD_y_B(3:end,:,:)-fD_y_B(1:end-2,:,:)).*aux_Bound(2:end-1,:,:);
            -fD_y_B(end-1,:,:).*aux_Bound(end,:,:)]+...
            aux_fD_z_B)...
            + (P.a_rp/(2*P.delta_x)).*...
            ([fU_x_B(:,2,:).*aux_Bound(:,1,:),...
            (fU_x_B(:,3:end,:)-fU_x_B(:,1:end-2,:)).*aux_Bound(:,2:end-1,:),...
            -fU_x_B(:,end-1,:).*aux_Bound(:,end,:)] + ...
            [fU_y_B(2,:,:).*aux_Bound(1,:,:);
            (fU_y_B(3:end,:,:)-fU_y_B(1:end-2,:,:)).*aux_Bound(2:end-1,:,:);
            -fU_y_B(end-1,:,:).*aux_Bound(end,:,:)]+...
            aux_fU_z_B);
        
    aux_vel_a = -(v_mm/(2*P.delta_x)).*...
            ([fD_x_a(:,2,:).*aux_Bound(:,1,:),...
            (fD_x_a(:,3:end,:)-fD_x_a(:,1:end-2,:)).*aux_Bound(:,2:end-1,:),...
            -fD_x_a(:,end-1,:).*aux_Bound(:,end,:)] + ...
            [fD_y_a(2,:,:).*aux_Bound(1,:,:);
            (fD_y_a(3:end,:,:)-fD_y_a(1:end-2,:,:)).*aux_Bound(2:end-1,:,:);
            -fD_y_a(end-1,:,:).*aux_Bound(end,:,:)]+...
            aux_fD_z_a);    
        
    aux_vel_c = -(v_mm/(2*P.delta_x)).*...
        ([fD_x_c(:,2,:).*aux_Bound(:,1,:),...
        (fD_x_c(:,3:end,:)-fD_x_c(:,1:end-2,:)).*aux_Bound(:,2:end-1,:),...
        -fD_x_c(:,end-1,:).*aux_Bound(:,end,:)] + ...
        [fD_y_c(2,:,:).*aux_Bound(1,:,:);
        (fD_y_c(3:end,:,:)-fD_y_c(1:end-2,:,:)).*aux_Bound(2:end-1,:,:);
        -fD_y_c(end-1,:,:).*aux_Bound(end,:,:)]+...
        aux_fD_z_c);
    
    
    dB_dt = aux_vel_B + P.psi_0*(f_nuc + f_sev + P.I_B*Omega+ J(1)*Omega2) - P.k_B*B;
    
    da_dt = aux_vel_a - P.k_A*a + P.I_A*Omega - f_nuc +J(2)*Omega2;
    
    dc_dt =  aux_vel_c - P.k_C*c + P.I_C*Omega - f_sev + J(3)*Omega2; 
        
    dvar_dt = [dB_dt(:);da_dt(:);dc_dt(:)];
end