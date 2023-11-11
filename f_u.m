function [f_u,f_d] = f_u(S,T_s,D,P)
    f_u = zeros(size(P.r));
    f_d = zeros(size(P.r));
    
    N = cross(S(T_s(:,2),:)-S(T_s(:,1),:),...
        S(T_s(:,3),:)-S(T_s(:,1),:));
    n = -N./sqrt(dot(N,N,2));
%     in the mesh the normal vectors point outside the spine, hence we take
%     the negative of the normal vector

    v = P.r-S(D(:,2),:);
    dist = D(:,1);
 

    f_u(D(:,4) == 1,:) = -dH(dist(D(:,4) == 1),P.alpha,P.beta)...
            .*v(D(:,4) == 1,:)./abs(dist(D(:,4) == 1));
    f_u(D(:,4) == 2,:) = -dH(dist(D(:,4) == 2),P.alpha,P.beta)...
            .*sign(dist(D(:,4) == 2)).*2.*n(D(D(:,4)==2,3),:);
    f_d(D(:,4) == 1,:) = v(D(:,4) == 1,:)./abs(dist(D(:,4) == 1));
    f_d(D(:,4) == 2,:) = 2.*sign(dist(D(:,4) == 2)).*n(D(D(:,4)==2,3),:);


end
                