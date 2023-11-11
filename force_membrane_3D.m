function f = force_membrane_3D(S,T_s,P)
    edges_s = edges(triangulation(double(T_s),S));
%     the spontaneous curvature is set to zero
    P.H_bar = 0;
    
    s1 = S(T_s(:,1),:);
    s2 = S(T_s(:,2),:);
    s3 = S(T_s(:,3),:);
    N = cross(s2-s1,s3-s1);
    N_l = sqrt(dot(N,N,2));
    n = N./N_l;
    A_all = N_l./2;
    %     bending contribution 
    f_b = zeros(size(S));
    
    
    H = zeros(size(S,1),1);
    
    for ind = 1:size(H,1)
        T1 = find(T_s(:,1) == ind);
        T2 = find(T_s(:,2) == ind);
        T3 = find(T_s(:,3) == ind);
        T_all = unique([T1;T2;T3]);
        A_i = sum(A_all(T_all))/3;
 
        e1 = find(edges_s(:,1) == ind);
        if ~isempty(e1)
            for l=1:length(e1)
                T_aux_i = [find(T_s(:,1) == edges_s(e1(l),1) & T_s(:,2) == edges_s(e1(l),2));
                    find(T_s(:,2) == edges_s(e1(l),1) & T_s(:,3) == edges_s(e1(l),2));
                    find(T_s(:,3) == edges_s(e1(l),1) & T_s(:,1) == edges_s(e1(l),2))];
                T_aux_j = [find(T_s(:,1) == edges_s(e1(l),2) & T_s(:,2) == edges_s(e1(l),1));
                    find(T_s(:,2) == edges_s(e1(l),2) & T_s(:,3) == edges_s(e1(l),1));
                    find(T_s(:,3) == edges_s(e1(l),2) & T_s(:,1) == edges_s(e1(l),1))];
                if ~isempty(T_aux_i) && ~isempty(T_aux_j)
                    l_ij = S(edges_s(e1(l),1),:) - S(edges_s(e1(l),2),:);
                    l_ij = sqrt(dot(l_ij,l_ij,2));
                    n_i = n(T_aux_i,:);
                    n_j = n(T_aux_j,:);
                    cos_theta = round(dot(n_i,n_j,2),10);
                    phi_ij = acos(cos_theta);
                    H(ind) = H(ind) + l_ij*phi_ij/(4*A_i);
                end
            end
        end

        e2 = find(edges_s(:,2) == ind);
        if ~isempty(e2)
            for l=1:length(e2)
                T_aux_i = [find(T_s(:,1) == edges_s(e2(l),2) & T_s(:,2) == edges_s(e2(l),1));
                    find(T_s(:,2) == edges_s(e2(l),2) & T_s(:,3) == edges_s(e2(l),1));
                    find(T_s(:,3) == edges_s(e2(l),2) & T_s(:,1) == edges_s(e2(l),1))];
                T_aux_j = [find(T_s(:,1) == edges_s(e2(l),1) & T_s(:,2) == edges_s(e2(l),2));
                    find(T_s(:,2) == edges_s(e2(l),1) & T_s(:,3) == edges_s(e2(l),2));
                    find(T_s(:,3) == edges_s(e2(l),1) & T_s(:,1) == edges_s(e2(l),2))];
                if ~isempty(T_aux_i) && ~isempty(T_aux_j)
                    l_ij = S(edges_s(e2(l),2),:) - S(edges_s(e2(l),1),:);
                    l_ij = sqrt(dot(l_ij,l_ij,2));
                    n_i = n(T_aux_i,:);
                    n_j = n(T_aux_j,:);
                    cos_theta = round(dot(n_i,n_j,2),10);
                    phi_ij = acos(cos_theta);
                    H(ind) = H(ind) + l_ij*phi_ij/(4*A_i);
                end
            end
        end
    end
    
    for ind = 1:size(S,1)
        e1 = find(edges_s(:,1) == ind);
        if ~isempty(e1)
            for l = 1:length(e1)
                T_aux_i1 = find(T_s(:,1) == edges_s(e1(l),1) & T_s(:,2) == edges_s(e1(l),2));
                T_aux_i2 = find(T_s(:,2) == edges_s(e1(l),1) & T_s(:,3) == edges_s(e1(l),2));
                T_aux_i3 = find(T_s(:,3) == edges_s(e1(l),1) & T_s(:,1) == edges_s(e1(l),2));
                T_aux_j1 = find(T_s(:,1) == edges_s(e1(l),2) & T_s(:,2) == edges_s(e1(l),1));
                T_aux_j2 = find(T_s(:,2) == edges_s(e1(l),2) & T_s(:,3) == edges_s(e1(l),1));
                T_aux_j3 = find(T_s(:,3) == edges_s(e1(l),2) & T_s(:,1) == edges_s(e1(l),1));
                v_ij = S(edges_s(e1(l),1),:) - S(edges_s(e1(l),2),:);
                l_ij = sqrt(dot(v_ij,v_ij,2));
                if l_ij == 0 
                    disp('error2')
                    break;
                end
                if ~isempty([T_aux_i1;T_aux_i2;T_aux_i3]) && ~isempty([T_aux_j1;T_aux_j2;T_aux_j3])
                    if ~isempty(T_aux_i1) 
                        n_i = n(T_aux_i1,:);
                        Si = S(T_s(T_aux_i1,1),:);
                        Sj = S(T_s(T_aux_i1,2),:);
                        Sk = S(T_s(T_aux_i1,3),:);
                        if ~isempty(T_aux_j1)  
                            n_j = n(T_aux_j1,:); 
                            Sl = S(T_s(T_aux_j1,3),:);
                        elseif ~isempty(T_aux_j2) 
                            n_j = n(T_aux_j2,:); 
                            Sl = S(T_s(T_aux_j2,1),:);
                        elseif ~isempty(T_aux_j3) 
                            n_j = n(T_aux_j3,:); 
                            Sl = S(T_s(T_aux_j3,2),:);
                        end
                    elseif ~isempty(T_aux_i2) 
                        n_i = n(T_aux_i2,:);
                        Si = S(T_s(T_aux_i2,2),:);
                        Sj = S(T_s(T_aux_i2,3),:);
                        Sk = S(T_s(T_aux_i2,1),:);
                        if ~isempty(T_aux_j1)  
                            n_j = n(T_aux_j1,:); 
                            Sl = S(T_s(T_aux_j1,3),:);
                        elseif ~isempty(T_aux_j2) 
                            n_j = n(T_aux_j2,:); 
                            Sl = S(T_s(T_aux_j2,1),:);
                        elseif ~isempty(T_aux_j3) 
                            n_j = n(T_aux_j3,:); 
                            Sl = S(T_s(T_aux_j3,2),:);
                        end
                    elseif ~isempty(T_aux_i3) 
                        n_i = n(T_aux_i3,:);
                        Si = S(T_s(T_aux_i3,3),:);
                        Sj = S(T_s(T_aux_i3,2),:);
                        Sk = S(T_s(T_aux_i3,1),:);
                        if ~isempty(T_aux_j1)  
                            n_j = n(T_aux_j1,:); 
                            Sl = S(T_s(T_aux_j1,3),:);
                        elseif ~isempty(T_aux_j2) 
                            n_j = n(T_aux_j2,:); 
                            Sl = S(T_s(T_aux_j2,1),:);
                        elseif ~isempty(T_aux_j3) 
                            n_j = n(T_aux_j3,:); 
                            Sl = S(T_s(T_aux_j3,2),:);
                        end
                    end
                    cos_theta = round(dot(n_i,n_j,2),10);
                    phi_ij = acos(cos_theta);
                    K_ij = phi_ij*v_ij./(2*l_ij);
                    v_jk = Sk - Sj;
                    w_jk = cross(Sk,Sj);
                    v_lj = Sj - Sl;
                    w_lj = cross(Sj,Sl);
                    H_ij = (cross(v_jk,w_jk) + cross(v_lj,w_lj))./4;
                    aux1 = cross(Si-Sj,Sk-Sj);
                    aux2 = cross(Si-Sj,Sl-Sj);

                    if sqrt(dot(aux1,aux1,2)) == 0 && sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_1 = [0 0 0];
                    elseif sqrt(dot(aux1,aux1,2)) == 0 
                        S_ij_1 = n_j*dot(Si-Sj,Sl-Sj,2)./sqrt(dot(aux2,aux2,2))./2;
                    elseif sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_1 =(n_i*dot(Si-Sj,Sk-Sj,2)./sqrt(dot(aux1,aux1,2)))./2;
                    else
                        S_ij_1 =  (n_i*dot(Si-Sj,Sk-Sj,2)./sqrt(dot(aux1,aux1,2)) + ...
                            n_j*dot(Si-Sj,Sl-Sj,2)./sqrt(dot(aux2,aux2,2)))./2;
                    end
              
                    aux1 = cross(Sj-Sk,Si-Sk);
                    aux2 = cross(Si-Sl,Sj-Sl);
                    if sqrt(dot(aux1,aux1,2)) == 0 && sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_2 = [0 0 0];
                    elseif sqrt(dot(aux1,aux1,2)) == 0 
                        S_ij_2 =  -(n_j*dot(Si-Sl,Sj-Sl,2)./sqrt(dot(aux2,aux2,2)));
                    elseif sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_2 =  -(n_i*dot(Sj-Sk,Si-Sk,2)./sqrt(dot(aux1,aux1,2)));
                    else
                        S_ij_2 =  -(n_i*dot(Sj-Sk,Si-Sk,2)./sqrt(dot(aux1,aux1,2)) + ...
                            n_j*dot(Si-Sl,Sj-Sl,2)./sqrt(dot(aux2,aux2,2)));
                    end
                    
                    f_b(ind,:) = f_b(ind,:) ...
                        -(P.kappa*(H(ind)-P.H_bar)+P.kappa*(H(edges_s(e1(l),2))-P.H_bar))*K_ij ...
                        +(P.kappa*(H(ind)-P.H_bar)*(H(ind)+P.H_bar)./3 + ...
                        2*P.kappa*(H(edges_s(e1(l),2))-P.H_bar)*(H(edges_s(e1(l),2))+P.H_bar)./3)*H_ij ...
                        -(P.kappa*(H(ind)-P.H_bar)*S_ij_1+P.kappa*(H(edges_s(e1(l),2))-P.H_bar)*S_ij_2);
                end
            end
        end
                          
        e1 = find(edges_s(:,2) == ind);
        if ~isempty(e1)
            for l = 1:length(e1)
                T_aux_i1 = find(T_s(:,1) == edges_s(e1(l),2) & T_s(:,2) == edges_s(e1(l),1));
                T_aux_i2 = find(T_s(:,2) == edges_s(e1(l),2) & T_s(:,3) == edges_s(e1(l),1));
                T_aux_i3 = find(T_s(:,3) == edges_s(e1(l),2) & T_s(:,1) == edges_s(e1(l),1));
                T_aux_j1 = find(T_s(:,1) == edges_s(e1(l),1) & T_s(:,2) == edges_s(e1(l),2));
                T_aux_j2 = find(T_s(:,2) == edges_s(e1(l),1) & T_s(:,3) == edges_s(e1(l),2));
                T_aux_j3 = find(T_s(:,3) == edges_s(e1(l),1) & T_s(:,1) == edges_s(e1(l),2));
                v_ij = S(edges_s(e1(l),2),:) - S(edges_s(e1(l),1),:);
                l_ij = sqrt(dot(v_ij,v_ij,2));
                if l_ij == 0 
                    disp('error1')
                    break;
                end
                
                if ~isempty([T_aux_i1;T_aux_i2;T_aux_i3]) && ~isempty([T_aux_j1;T_aux_j2;T_aux_j3])
                    if ~isempty(T_aux_i1) 
                        n_i = n(T_aux_i1,:);
                        Si = S(T_s(T_aux_i1,1),:);
                        Sj = S(T_s(T_aux_i1,2),:);
                        Sk = S(T_s(T_aux_i1,3),:);
                        if ~isempty(T_aux_j1)  
                            n_j = n(T_aux_j1,:); 
                            Sl = S(T_s(T_aux_j1,3),:);
                        elseif ~isempty(T_aux_j2) 
                            n_j = n(T_aux_j2,:); 
                            Sl = S(T_s(T_aux_j2,1),:);
                        elseif ~isempty(T_aux_j3) 
                            n_j = n(T_aux_j3,:); 
                            Sl = S(T_s(T_aux_j3,2),:);
                        end
                    elseif ~isempty(T_aux_i2) 
                        n_i = n(T_aux_i2,:);
                        Si = S(T_s(T_aux_i2,2),:);
                        Sj = S(T_s(T_aux_i2,3),:);
                        Sk = S(T_s(T_aux_i2,1),:);
                        if ~isempty(T_aux_j1)  
                            n_j = n(T_aux_j1,:); 
                            Sl = S(T_s(T_aux_j1,3),:);
                        elseif ~isempty(T_aux_j2) 
                            n_j = n(T_aux_j2,:); 
                            Sl = S(T_s(T_aux_j2,1),:);
                        elseif ~isempty(T_aux_j3) 
                            n_j = n(T_aux_j3,:); 
                            Sl = S(T_s(T_aux_j3,2),:);
                        end
                    elseif ~isempty(T_aux_i3) 
                        n_i = n(T_aux_i3,:);
                        Si = S(T_s(T_aux_i3,3),:);
                        Sj = S(T_s(T_aux_i3,2),:);
                        Sk = S(T_s(T_aux_i3,1),:);
                        if ~isempty(T_aux_j1)  
                            n_j = n(T_aux_j1,:); 
                            Sl = S(T_s(T_aux_j1,3),:);
                        elseif ~isempty(T_aux_j2) 
                            n_j = n(T_aux_j2,:); 
                            Sl = S(T_s(T_aux_j2,1),:);
                        elseif ~isempty(T_aux_j3) 
                            n_j = n(T_aux_j3,:); 
                            Sl = S(T_s(T_aux_j3,2),:);
                        end
                    end
                    cos_theta = round(dot(n_i,n_j,2),10);
                    phi_ij = acos(cos_theta);
                    K_ij = phi_ij*v_ij./(2*l_ij);
                    v_jk = Sk - Sj;
                    w_jk = cross(Sk,Sj);
                    v_lj = Sj - Sl;
                    w_lj = cross(Sj,Sl);
                    H_ij = (cross(v_jk,w_jk) + cross(v_lj,w_lj))./4;

                    aux1 = cross(Si-Sj,Sk-Sj);
                    aux2 = cross(Si-Sj,Sl-Sj);
                    if sqrt(dot(aux1,aux1,2)) == 0 && sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_1 = [0 0 0];
                    elseif sqrt(dot(aux1,aux1,2)) == 0 
                        S_ij_1 = n_j*dot(Si-Sj,Sl-Sj,2)./sqrt(dot(aux2,aux2,2))./2;
                    elseif sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_1 =  (n_i*dot(Si-Sj,Sk-Sj,2)./sqrt(dot(aux1,aux1,2)))./2;
                    else
                        S_ij_1 =  (n_i*dot(Si-Sj,Sk-Sj,2)./sqrt(dot(aux1,aux1,2)) + ...
                            n_j*dot(Si-Sj,Sl-Sj,2)./sqrt(dot(aux2,aux2,2)))./2;
                    end
                    
                    aux1 = cross(Sj-Sk,Si-Sk);
                    aux2 = cross(Si-Sl,Sj-Sl);
                    if sqrt(dot(aux1,aux1,2)) == 0 && sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_2 = [0 0 0];
                    elseif sqrt(dot(aux1,aux1,2)) == 0 
                        S_ij_2 =  -(n_j*dot(Si-Sl,Sj-Sl,2)./sqrt(dot(aux2,aux2,2)));
                    elseif sqrt(dot(aux2,aux2,2)) == 0  
                        S_ij_2 =  -(n_i*dot(Sj-Sk,Si-Sk,2)./sqrt(dot(aux1,aux1,2)));
                    else
                        S_ij_2 =  -(n_i*dot(Sj-Sk,Si-Sk,2)./sqrt(dot(aux1,aux1,2)) + ...
                            n_j*dot(Si-Sl,Sj-Sl,2)./sqrt(dot(aux2,aux2,2)));
                    end
                    
                    
                    f_b(ind,:) = f_b(ind,:) ...
                        -(P.kappa*(H(ind)-P.H_bar)+P.kappa*(H(edges_s(e1(l),1))-P.H_bar))*K_ij ...
                        +(P.kappa*(H(ind)-P.H_bar)*(H(ind)+P.H_bar)./3 + ...
                        2*P.kappa*(H(edges_s(e1(l),1))-P.H_bar)*(H(edges_s(e1(l),1))+P.H_bar)./3)*H_ij ...
                        -(P.kappa*(H(ind)-P.H_bar)*S_ij_1+P.kappa*(H(edges_s(e1(l),1))-P.H_bar)*S_ij_2); 
                end
            end
        end
    end
        
    f = f_b;
end
 
