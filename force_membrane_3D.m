function f = force_membrane_3D(S,T_s,P)
    
%     bending engergy contribution
    f_b = zeros(size(S));
       
    edges_s = edges(triangulation(double(T_s),S));
    for k = 1:size(edges_s,1)
        aux1 = edges_s(k,1);
        aux2 = edges_s(k,2);
        si = [];
        sip1 = [];
        sip2 = [];
        sj = [];
        sjp1 = [];
        sjp2 = [];
        
        aux_Ti = find(T_s(:,1)==aux1);
        if ~isempty(aux_Ti)
            aux_Tip1 = find(T_s(aux_Ti,2)==aux2);
            if ~isempty(aux_Tip1)
                si = T_s(aux_Ti(aux_Tip1),1);
                sip1 = T_s(aux_Ti(aux_Tip1),2);
                sip2 = T_s(aux_Ti(aux_Tip1),3);
                sj = sip1;
                sjp1 = si;
                
  %                 case s_i = s_j
                aux_j = find(T_s(:,1)==aux2);
                if ~isempty(aux_j)
                    aux_jp = find(T_s(aux_j,2)==aux1);
                    if ~isempty(aux_jp)
                        sjp2 = T_s(aux_j(aux_jp),3);
                    end
                else
  %                 case s_i = s_{j+1}
                    aux_j = find(T_s(:,2)==aux2);
                    if ~isempty(aux_j)
                        aux_jp = find(T_s(aux_j,3)==aux1);
                        if ~isempty(aux_jp)
                            sjp2 = T_s(aux_j(aux_jp),1);
                        end
                    else
         %                 case s_i = s_{j+2}
                        aux_j = find(T_s(:,3)==aux2);
                        if ~isempty(aux_j)
                            aux_jp = find(T_s(aux_j,1)==aux1);
                            if ~isempty(aux_jp)
                                sjp2 = T_s(aux_j(aux_jp),2);
                            end
                        end
                    end
                end
            else
                aux_Ti = find(T_s(:,2)==aux1);
                if ~isempty(aux_Ti)
                    aux_Tip1 = find(T_s(aux_Ti,3)==aux2);
                    if ~isempty(aux_Tip1)
                        si = T_s(aux_Ti(aux_Tip1),2);
                        sip1 = T_s(aux_Ti(aux_Tip1),3);
                        sip2 = T_s(aux_Ti(aux_Tip1),1);
                        sj = sip1;
                        sjp1 = si;
                        %                 case s_{i+1} = s_{j}
                        aux_j = find(T_s(:,1)==aux2);
                        if ~isempty(aux_j)
                            aux_jp = find(T_s(aux_j,2)==aux1);
                            if ~isempty(aux_jp)
                                sjp2 = T_s(aux_j(aux_jp),3);
                            end
                        else
                            %                 case s_{i+1} = s_{j+1}
                            aux_j = find(T_s(:,2)==aux2);
                            if ~isempty(aux_j)
                                aux_jp = find(T_s(aux_j,3)==aux1);
                                if ~isempty(aux_jp)
                                    sjp2 = T_s(aux_j(aux_jp),1);
                                end
                            else
                                %                 case s_{i+1} = s_{j+2}
                                aux_j = find(T_s(:,3)==aux2);
                                if ~isempty(aux_j)
                                    aux_jp = find(T_s(aux_j,1)==aux1);
                                    if ~isempty(aux_jp)
                                        sjp2 = T_s(aux_j(aux_jp),2);
                                    end
                                end
                            end
                        end
                    else
                        aux_Ti = find(T_s(:,3)==aux1);
                        if ~isempty(aux_Ti)
                            aux_Tip1 = find(T_s(aux_Ti,1)==aux2);
                            if ~isempty(aux_Tip1)
                                si = T_s(aux_Ti(aux_Tip1),3);
                                sip1 = T_s(aux_Ti(aux_Tip1),1);
                                sip2 = T_s(aux_Ti(aux_Tip1),2);
                                sj = sip1;
                                sjp1 = si;
            %                 case s_{i+2} = s_{j}
                                aux_j = find(T_s(:,1)==aux2);
                                if ~isempty(aux_j)
                                    aux_jp = find(T_s(aux_j,2)==aux1);
                                    if ~isempty(aux_jp)
                                        sjp2 = T_s(aux_j(aux_jp),3);
                                    end
                                else
                                    %                 case s_{i+2} = s_{j+1}
                                    aux_j = find(T_s(:,2)==aux2);
                                    if ~isempty(aux_j)
                                        aux_jp = find(T_s(aux_j,3)==aux1);
                                        if ~isempty(aux_jp)
                                            sjp2 = T_s(aux_j(aux_jp),1);
                                        end
                                    else
                                        %                 case s_{i+2} = s_{j+2}
                                        aux_j = find(T_s(:,3)==aux2);
                                        if ~isempty(aux_j)
                                            aux_jp = find(T_s(aux_j,1)==aux1);
                                            if ~isempty(aux_jp)
                                                sjp2 = T_s(aux_j(aux_jp),2);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
%         changing i for j and checking all the cases 
        if isempty(si)
            aux1 = edges_s(k,2);
            aux2 = edges_s(k,1);
            aux_Tj = find(T_s(:,1)==aux1);
            if ~isempty(aux_Tj)
                aux_Tjp1 = find(T_s(aux_Tj,2)==aux2);
                if ~isempty(aux_Tjp1)
                    sj = T_s(aux_Tj(aux_Tjp1),1);
                    sjp1 = T_s(aux_Tj(aux_Tjp1),2);
                    sjp2 = T_s(aux_Tj(aux_Tjp1),3);
                    si = sjp1;
                    sip1 = sj;

                    aux_i = find(T_s(:,1)==aux2);
                    if ~isempty(aux_i)
                        aux_ip = find(T_s(aux_i,2)==aux1);
                        if ~isempty(aux_ip)
                            sip2 = T_s(aux_i(aux_ip),3);
                        end
                    else
                        aux_i = find(T_s(:,2)==aux2);
                        if ~isempty(aux_i)
                            aux_ip = find(T_s(aux_i,3)==aux1);
                            if ~isempty(aux_ip)
                                sip2 = T_s(aux_i(aux_ip),1);
                            end
                        else
                            aux_i = find(T_s(:,3)==aux2);
                            if ~isempty(aux_i)
                                aux_ip = find(T_s(aux_i,1)==aux1);
                                if ~isempty(aux_ip)
                                    sip2 = T_s(aux_i(aux_ip),2);
                                end
                            end
                        end
                    end
                else
                    aux_Tj = find(T_s(:,2)==aux1);
                    if ~isempty(aux_Tj)
                        aux_Tjp1 = find(T_s(aux_Tj,3)==aux2);
                        if ~isempty(aux_Tjp1)
                            sj = T_s(aux_Tj(aux_Tjp1),1);
                            sjp1 = T_s(aux_Tj(aux_Tjp1),2);
                            sjp2 = T_s(aux_Tj(aux_Tjp1),3);
                            si = sjp1;
                            sip1 = sj;

                            aux_i = find(T_s(:,1)==aux2);
                            if ~isempty(aux_i)
                                aux_ip = find(T_s(aux_i,2)==aux1);
                                if ~isempty(aux_ip)
                                    sip2 = T_s(aux_i(aux_ip),3);
                                end
                            else
                                aux_i = find(T_s(:,2)==aux2);
                                if ~isempty(aux_i)
                                    aux_ip = find(T_s(aux_i,3)==aux1);
                                    if ~isempty(aux_ip)
                                        sip2 = T_s(aux_i(aux_ip),1);
                                    end
                                else
                                    aux_i = find(T_s(:,3)==aux2);
                                    if ~isempty(aux_i)
                                        aux_ip = find(T_s(aux_i,1)==aux1);
                                        if ~isempty(aux_ip)
                                            sip2 = T_s(aux_i(aux_ip),2);
                                        end
                                    end
                                end
                            end
                        else
                            aux_Tj = find(T_s(:,3)==aux1);
                            if ~isempty(aux_Tj)
                                aux_Tjp1 = find(T_s(aux_Tj,1)==aux2);
                                if ~isempty(aux_Tjp1)
                                    sj = T_s(aux_Tj(aux_Tjp1),1);
                                    sjp1 = T_s(aux_Tj(aux_Tjp1),2);
                                    sjp2 = T_s(aux_Tj(aux_Tjp1),3);
                                    si = sjp1;
                                    sip1 = sj;

                                    aux_i = find(T_s(:,1)==aux2);
                                    if ~isempty(aux_i)
                                        aux_ip = find(T_s(aux_i,2)==aux1);
                                        if ~isempty(aux_ip)
                                            sip2 = T_s(aux_i(aux_ip),3);
                                        end
                                    else
                                        aux_i = find(T_s(:,2)==aux2);
                                        if ~isempty(aux_i)
                                            aux_ip = find(T_s(aux_i,3)==aux1);
                                            if ~isempty(aux_ip)
                                                sip2 = T_s(aux_i(aux_ip),1);
                                            end
                                        else
                                            aux_i = find(T_s(:,3)==aux2);
                                            if ~isempty(aux_i)
                                                aux_ip = find(T_s(aux_i,1)==aux1);
                                                if ~isempty(aux_ip)
                                                    sip2 = T_s(aux_i(aux_ip),2);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

        end
%         if the edge is shared by two triangles in the mesh, calculate the
%         force generated by the bending force in each vertex of both
%         triangles
        if  ~isempty(sjp2) && ~isempty(sip2)
            Ni = cross(S(sip1,:)-S(si,:),S(sip2,:)-S(si,:));
            length_Ni = sqrt(dot(Ni,Ni,2));
            ni = Ni./length_Ni;
            Nj = cross(S(sjp1,:)-S(sj,:),S(sjp2,:)-S(sj,:));
            length_Nj = sqrt(dot(Nj,Nj,2));
            nj = Nj./length_Nj;
            aux_n = dot(ni,nj,2);
            
            f_b(si,:) = f_b(si,:) + cross(nj-aux_n.*ni,S(sip2,:)-S(sip1,:))./length_Ni ...
                + cross(ni-aux_n.*nj,S(sj,:)-S(sjp2,:))./length_Nj;
            
            f_b(sip1,:) = f_b(sip1,:) + cross(nj-aux_n.*ni,S(si,:)-S(sip2,:))./length_Ni ...
                + cross(ni-aux_n.*nj,S(sjp2,:)-S(sjp1,:))./length_Nj;
            
            f_b(sip2,:) = f_b(sip2,:) + cross(nj-aux_n.*ni,S(sip1,:)-S(si,:))./length_Ni;
            
            f_b(sjp2,:) = f_b(sjp2,:) + cross(ni-aux_n.*nj,S(sjp1,:)-S(sj,:))./length_Nj;
        end 
        
    end
    f =  P.kappa*f_b;
end
