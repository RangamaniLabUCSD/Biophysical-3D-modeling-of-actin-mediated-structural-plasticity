function d = dist_r_s_3D(S,T_s,P)
  
    d = zeros(P.n_x3,4);%[distance index_S index_T (S==1,T==2)]
    N = cross(S(T_s(:,2),:) - S(T_s(:,1),:),S(T_s(:,3),:) - S(T_s(:,1),:));%normal vector
    n = N./sqrt(dot(N,N,2));
    for k = 1:P.n_x3
%         calculate distance to a vertex
        v = P.r(k,:) - S;
        d_v = sqrt(dot(v,v,2));
        [min_val_S,ind_min_S] = min(d_v);
        aux1 = find(T_s(:,1) == ind_min_S);
        aux2 = find(T_s(:,2) == ind_min_S);
        aux3 = find(T_s(:,3) == ind_min_S);
        aux_all = [aux1;aux2;aux3];
        
        
%         calculate distance to a the plane spanned by a triangular face of
%         the mesh
        aux_d = abs(dot(P.r(k,:)- S(T_s(aux_all,1),:),n(aux_all,:),2));
        [min_val_p,ind_min_p] = min(aux_d);
                
        if min_val_S <= min_val_p
            d(k,:) = [min_val_S ind_min_S aux_all(ind_min_p) 1];
        else
            d(k,:) = [min_val_p ind_min_S aux_all(ind_min_p) 2];
        end
%         check sign function 
        aux_v = P.r(k,:) - S(T_s(aux_all(ind_min_p),1),:);
        aux_dir = acos(dot(aux_v,n(aux_all(ind_min_p),:),2)./sqrt(dot(aux_v,aux_v,2)));
        if aux_dir > pi/2
            d(k,1) = -d(k,1);
        end
    end
end

     
        