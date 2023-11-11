function f = force_actin_3D(S,T_s,D,B,alpha,beta,P)
    f = zeros(size(S));
    
    N = cross(S(T_s(:,2),:)-S(T_s(:,1),:),...
        S(T_s(:,3),:)-S(T_s(:,1),:));
    N_length = sqrt(dot(N,N,2));
    n = -N./N_length;
    ds1 = S(T_s(:,3),:)-S(T_s(:,2),:); 
    ds2 = S(T_s(:,1),:)-S(T_s(:,3),:); 
    ds3 = S(T_s(:,2),:)-S(T_s(:,1),:); 
    
    for k = 1:size(S,1)
        v =P.r - S(k,:);
%        for distance to a vertex 
        v_index = find(D(:,4) == 1 & D(:,2) == k);
%     for distance to a triangular face    
        T_index_s1 = find(T_s(:,1) == k);
        T_index_s2 = find(T_s(:,2) == k);
        T_index_s3 = find(T_s(:,3) == k);
        
        r_T_index_s1 = [];
        r_T_index_s2 = [];
        r_T_index_s3 = [];
        if ~isempty(T_index_s1)
            for l = 1:length(T_index_s1)
                aux = find(D(:,3)==T_index_s1(l) & D(:,4) == 2);
                r_T_index_s1 = [r_T_index_s1;
                    aux T_index_s1(l)*ones(size(aux))];
            end
        end
        if ~isempty(T_index_s2)
            for l = 1:length(T_index_s2)
                aux = find(D(:,3)==T_index_s2(l) & D(:,4) == 2);
                r_T_index_s2 = [r_T_index_s2;
                    aux T_index_s2(l)*ones(size(aux))];
            end
        end
        if ~isempty(T_index_s3)
            for l = 1:length(T_index_s3)
                aux = find(D(:,3)==T_index_s3(l) & D(:,4) == 2);
                r_T_index_s3 = [r_T_index_s3;
                    aux T_index_s3(l)*ones(size(aux))];
            end
        end
%  get the contribution of the force to the vertex       
        aux_v = zeros(1,3);
        if ~isempty(v_index)

           aux_v = (P.delta_x^3).*sum(dH(D(v_index,1),alpha,beta).*B(v_index).*...
                (-v(v_index,:))./abs(D(v_index,1)));
        end
        
        aux_T1 = zeros(1,3);
        if ~isempty(r_T_index_s1)

            aux_T1  = (P.delta_x^3).*sum(dH(D(r_T_index_s1(:,1),1),alpha,beta).*...
                B(r_T_index_s1(:,1)).*...
                2.*sign(D(r_T_index_s1(:,1),1)).*(...
                -n(r_T_index_s1(:,2),:) + cross(...
                v(r_T_index_s1(:,1),:) - ...
                dot(v(r_T_index_s1(:,1),:),n(r_T_index_s1(:,2),:),2).*n(r_T_index_s1(:,2)),...
                ds1(r_T_index_s1(:,2),:))./N_length(r_T_index_s1(:,2))));
        end
            
        aux_T2 = zeros(1,3);
        if ~isempty(r_T_index_s2)

        aux_T2 =  (P.delta_x^3).*sum(dH(D(r_T_index_s2(:,1),1),alpha,beta).*...
                B(r_T_index_s2(:,1)).*...
                2.*sign(D(r_T_index_s2(:,1),1)).*...
                cross(v(r_T_index_s2(:,1),:)-...
                dot(n(r_T_index_s2(:,2),:),v(r_T_index_s2(:,1),:),2).*n(r_T_index_s2(:,2),:),...
                ds2(r_T_index_s2(:,2),:))./N_length(r_T_index_s2(:,2),:));
        end
        
        aux_T3 = zeros(1,3);
        if ~isempty(r_T_index_s3)

            aux_T3 =  (P.delta_x^3).*sum(dH(D(r_T_index_s3(:,1),1),alpha,beta).*...
                B(r_T_index_s3(:,1)).*...
                2.*sign(D(r_T_index_s3(:,1),1)).*...
                cross(v(r_T_index_s3(:,1),:)-...
                dot(n(r_T_index_s3(:,2),:),v(r_T_index_s3(:,1),:),2).*n(r_T_index_s3(:,2),:),...
                ds3(r_T_index_s3(:,2),:))./N_length(r_T_index_s3(:,2),:));
        end
        
        f(k,:) = -(aux_T1 + aux_T2 + aux_T3 +aux_v);%
    end
end
        
  
        