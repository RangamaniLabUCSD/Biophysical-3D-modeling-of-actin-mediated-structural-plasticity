close all
clear 

load('Initial_spine.mat')


P.index_fix = [find(S(:,1) >= P.s_th);
            find(S(:,2) <= -P.s_th);
            find(S(:,1) <= -P.s_th);
            find(S(:,2) >= P.s_th)];
P.index_fix = intersect(P.index_fix,find(S(:,3) <= P.neck_ini/2));

ll_d_aux = [0.280236090138034,0.00557530226794533,0.881620479908731];
dist = sqrt(dot(S - ll_d_aux,S - ll_d_aux,2));
dist_ind = find(dist == min(dist));
ll_d = S(dist_ind,:);

% initiate time
P.t_ini = 0;%s, initial time 
P.t_end = 5*60;%s, final time 
P.t_si = P.t_ini;%s, initial stimulation time 
P.t_se = P.t_si+60;%s, final time of LTP influx
aux_t = P.t_ini:P.delta_t:P.t_end;

% variables
save_aux = 10;
S_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);%shape
T_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);%triangulation
var_save = zeros(P.t_end/(save_aux*P.delta_t)+1,3*P.n_x3);%[B,A,C]
v_mp = zeros(P.t_end/(save_aux*P.delta_t)+1,1);%protrusion velocity

l=1;
S_save{l,1} = S;
T_save{l,1} = T_s;
var_save(l,:) = var_0;
v_mp(l) = v_mm;

for k = 2:length(aux_t)
    if aux_t(k) > P.t_si && aux_t(k)<=P.t_se
        [S,var_0,ll_d,v_mm] = solve_system_threshold_3D(S,T_s,var_0,aux_t(k-1),P.J,ll_d,v_mm,P); 
        [T_s,S] = remeshing(int32(T_s),S,int32(P.index_fix),P.delta_S,int32(3));
        P.index_fix = [find(S(:,1) >= P.s_th);
            find(S(:,2) <= -P.s_th);
            find(S(:,1) <= -P.s_th);
            find(S(:,2) >= P.s_th)];
        P.index_fix = intersect(P.index_fix,find(S(:,3) <= P.neck_ini/2));
        S(P.index_fix,3) = 0; 
    else
        [S,var_0,ll_d,v_mm] = solve_system_threshold_3D(S,T_s,var_0,aux_t(k-1),P.J0,ll_d,v_mm,P); 
        [T_s,S] = remeshing(int32(T_s),S,int32(P.index_fix),P.delta_S,int32(3));
        P.index_fix = [find(S(:,1) >= P.s_th);
            find(S(:,2) <= -P.s_th);
            find(S(:,1) <= -P.s_th);
            find(S(:,2) >= P.s_th)];
        P.index_fix = intersect(P.index_fix,find(S(:,3) <= P.neck_ini/2));
        S(P.index_fix,3) = 0; 
    end
    if mod(k,save_aux) == 1
        
        disp(aux_t(k))
        l=l+1;
        S_save{l,1} = S;
        T_save{l,1} = T_s;
        var_save(l,:) = var_0;
        v_mp(l) = v_mm;
       
    end
   
    
end      

figure
% plot barbed ends
save_aux = 10;
color_p = pink;
color_p = color_p(end:-1:1,:);
l = fix(k/save_aux);
P.x = min(P.X(:)):P.delta_x:max(P.X(:));
xx = find(P.x < 0,1,'last');
XX = zeros(length(P.x));
YY = XX;
ZZ = XX;
XX(:,:)= P.X(:,xx,:);
YY(:,:)= P.Y(:,xx,:);
ZZ(:,:)= P.Z(:,xx,:);
set(gcf, 'Position',  [50, 50, 400, 400])
h1 = axes('Position', [0.1 0.26 0.8 0.7]);
aux_p = reshape(var_save(l,1:P.n_x3),P.n_x,P.n_x,P.n_x);
aux_pp = zeros(P.n_x);
aux_pp(:,:) = aux_p(xx,:,:);
scatter3(XX(:),YY(:),ZZ(:),50,aux_pp(:),'filled','s')
hold on 
trimesh(T_save{l,1},S_save{l,1}(:,1),S_save{l,1}(:,2),S_save{l,1}(:,3),'edgecolor','k','linewidth',1.5)
trimesh(T_save{1,1},S_save{1,1}(:,1),S_save{1,1}(:,2),S_save{1,1}(:,3),'edgecolor',[0.5 0.5 0.5],'linewidth',1)
view(90,0)
h = colorbar;
set(get(h,'label'),'string','[# Barbed ends/{\mu}m^3]');
colormap(h1,color_p)
h.Position = [0.15 0.12 0.7 0.02];
h.Orientation = 'horizontal';
caxis([0 15])
axis('square')
set(gcf,'color','w')
ylabel('[{\mu}m]')
zlabel('[{\mu}m]')
set(gca,'FontSize',14)
axis([-P.delta_S P.delta_S min(P.Y(:)) max(P.Y(:)) min(P.Z(:)) max(P.Z(:))])

