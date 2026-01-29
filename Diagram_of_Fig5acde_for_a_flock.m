function Diagram_of_Fig5acde_for_a_flock
addpath(genpath('BiMat'))

load group_05.mat
tracks_filt = tracks;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate Frame_matrix, tracks_filt for a flock
BirdIDs=unique(tracks(:,1));
T=unique(tracks(:,5));
Frame_matrix = zeros(length(BirdIDs),length(T));
for i = 1 : length(BirdIDs)
    Frame_matrix(i,:) = find(tracks(:,1)==BirdIDs(i));
    if sum(tracks(Frame_matrix(i,:),5)-T)~=0
       error('Error') 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute pairwise leader-follower relationship matrix

findRou = 'max';
tau_threshold = 0.1;
correlation_threshold = -1;
[Delay_net,One_LFnet,~,~] = Mapping_Leader_follow_network_anis_factor_simplified(Frame_matrix,tracks,findRou,tau_threshold,correlation_threshold);

[Rank,turn_interval]=GetRank(tracks_filt(Frame_matrix(:),:),tau_threshold,correlation_threshold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-time & 2-time modularity analysis

anis_factor = 0;
numBird = 70;

[all_module_bird,all_module_net] = Generate_Modularity(One_LFnet);
[all_module_bird_2times,module_bird_2times] = Generate_2times_Modules(One_LFnet,all_module_bird);

[~,Cr]=global_reaching_centrality(sign(abs(all_module_net))');
[~,sort_Cr] = sort(Cr,'descend');

load('roma.mat');
roma = roma(ceil(linspace(1,size(roma,1),size(all_module_net,1))),:);

tag_m1times = zeros(1,length(Rank));
for i = 1 : length(all_module_bird)
    tag_m1times(all_module_bird{i}) = i;
end
color_m1times = roma(tag_m1times,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot Fig.5cde
all_MS_type = {'max','min','middle','random'};
next_tau_slice = 50;
is_retina = 0;
consider_ego_motion = 'yes';
forward_or_omni = 'forward';

[all_MS_net,aveKij_Tree,aveKij_Tree_t,S_pass2Node,corr_S_T,ave_Cr,corr_Cr_T,std_tag_of_diffRoot,std_tag_of_diffRoot_t,std_spatial_of_diffRoot,std_spatial_of_diffRoot_t] = for_a_track( ...
    forward_or_omni,consider_ego_motion,is_retina,next_tau_slice,all_MS_type,Frame_matrix,tracks_filt,turn_interval,all_module_bird,all_module_bird_2times);

colors = [hex2rgb('0072BD');hex2rgb('D95319');hex2rgb('EDB120');hex2rgb('7E2F8E')];

figure;
set(gcf,'Position',[145 617 830 220])

subplot('Position',[0.07    0.14    0.25    0.78])
yline(0.5,'--')
X = aveKij_Tree_t';
box on;
violins = violinplot(X(:),[repmat({'1max-'},size(X,1),1);repmat({'2min-'},size(X,1),1); ...
    repmat({'3mid-'},size(X,1),1);repmat({'4rand-'},size(X,1),1)], ...
    'PointSize',6,'ShowMean',false);
for i  = 1 : 4
    violins(i).ViolinColor = colors(i,:);
    violins(i).PointEdgeColor = 'n';
    violins(i).PointColor = colors(i,:);
end
set(gca,'fontsize',14,'XTickLabelRotation',0,'YTick',[0.2 0.4 0.5 0.6 0.8])
ylim([0.2 0.8])
xlim([0.6 4.4])
ylabel('$k^{\rm{rank}}_{ij}(t,\tau_f)$','Interpreter','latex')


rows = 2;
cols = 2;
sub_pos = subplot_in_given_position([0.40 0.14 0.25 0.78],rows,cols);

% subplot('Position',[0.40    0.14    0.25    0.78]);hold on;
for i = 1 : length(all_MS_type)
    subplot('Position',sub_pos(i,:));hold on;
    x = turn_interval;
    y = ave_Cr(i,:);
    scatter(x,y,20,'MarkerEdgeColor','none','MarkerFaceColor',colors(i,:))
    box on;
    % [BS,xs]=func_cal_rlowess_bootstrap(x,y,10);
    % mean_BS=mean(BS,1);
    % Prctile_BS = prctile(BS,[3  97],1)';
    % h_area=area(xs',[Prctile_BS(:,1), (Prctile_BS(:,2)-Prctile_BS(:,1)) ],'FaceAlpha',0.5);
    % h_area(1).FaceColor = 'none';
    % h_area(2).FaceColor = hex2rgb('F7BECB');
    % h_area(1).EdgeColor = 'none';
    % h_area(2).EdgeColor = 'none';
    % h = plot(xs,mean_BS','-','linewidth',2,'color','r');

    xlim([0 3])
    ylim([0 0.1])
    grid on
    set(gca,'FontSize',14)
    if i==1 || i == 2
        set(gca,'XTickLabel','')
    end
    if i==2 || i == 4
        set(gca,'YTickLabel','')
    end
    if i==1 || i == 3
        ylabel('$\left\langle C_r^i\right\rangle_t$','Interpreter','latex')
    end
    if i==3 || i == 4
        xlabel('{\it T_i} (s)')
    end
end

% subplot('Position',[0.74    0.14    0.25    0.78])
rows = 2;
cols = 1;
sub_pos = subplot_in_given_position([0.74    0.14    0.25    0.78],rows,cols);

subplot('Position',sub_pos(1,:));
X = std_tag_of_diffRoot_t';
box on;
violins = violinplot(X(:),[repmat({'1max-'},size(X,1),1);repmat({'2min-'},size(X,1),1); ...
    repmat({'3mid-'},size(X,1),1);repmat({'4rand-'},size(X,1),1)], ...
    'PointSize',6,'ShowMean',false);
for i  = 1 : 4
    violins(i).ViolinColor = colors(i,:);
    violins(i).PointEdgeColor = 'n';
    violins(i).PointColor = colors(i,:);
end
set(gca,'fontsize',14,'XTickLabelRotation',0,'XTickLabel','','YTick',[0:0.1:0.3])
xlim([0.6 4.4])
ylim([0 0.3])
ylabel('$\sigma^{\rm{M}}_{\rm{ST}}(t)$','Interpreter','latex')

subplot('Position',sub_pos(2,:));
X = std_spatial_of_diffRoot_t';
box on;
violins = violinplot(X(:),[repmat({'1max-'},size(X,1),1);repmat({'2min-'},size(X,1),1); ...
    repmat({'3mid-'},size(X,1),1);repmat({'4rand-'},size(X,1),1)], ...
    'PointSize',6,'ShowMean',false);
for i  = 1 : 4
    violins(i).ViolinColor = colors(i,:);
    violins(i).PointEdgeColor = 'n';
    violins(i).PointColor = colors(i,:);
end
set(gca,'fontsize',14,'XTickLabelRotation',0,'YTick',[0:0.1:0.3])
xlim([0.6 4.4])
ylim([0 0.3])
ylabel('$\sigma^{\rm{S}}_{\rm{ST}}(t)$','Interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot Fig.5a

forward_or_omni = 'forward';
consider_ego_motion = 'yes';
is_retina = 0;
[all_MS_net_max,ave_Gr_maxMS,Gr_t_maxMS] = Cal_1_MSnet_for_a_track(forward_or_omni,consider_ego_motion,is_retina,'max',Frame_matrix,tracks_filt);


fig = figure;
set(gcf,'Position',[ 9    61   750   450])

rows = 1;
cols = 2;
sub_pos = subplot_in_given_position([0.02 0.36 0.6 0.6],rows,cols);
frame_t = [30];
BD=ceil([min(tracks(:,2))-1 max(tracks(:,2))+1 ...
    min(tracks(:,3))-1 max(tracks(:,3))+1 ...
    min(tracks(:,4))-1 max(tracks(:,4))+1]);
numBird = size(Frame_matrix,1);

red_color = [hex2rgb('FFCCCC');hex2rgb('F08080');hex2rgb('FF0000');hex2rgb('8B0000')];
red_color = turbo(4);
red_color(2,:) = [255 0 0]/255;
red_color(3,:) = hex2rgb('ffa500');
red_color(4,:) = hex2rgb('0072BD');



subplot('Position',sub_pos(2,:))

MS_net = squeeze(all_MS_net_max(frame_t(1),:,:));
[Gr_t,Cr_t,dis_t] = global_reaching_centrality(sign(MS_net));
g = digraph(MS_net);
[bins,binsize]= conncomp(g,'Type','weak');
root_node = zeros(1,length(binsize));
for i = 1 : length(binsize)
    index = find(bins==i);
    [~,index1] = max(Cr_t(index));
    root_node(i) = index(index1);

    TR{i} = shortestpathtree(g,root_node(i));
    rNodes = TR{i}.Edges.EndNodes(:,2)';

    temp1 = rNodes(find(MS_net(rNodes,root_node(i))))';
    TR{i} = addedge(TR{i},temp1,repmat(root_node(i),length(temp1),1),ones(length(temp1),1));

    mean_tag_t(i) = mean(tag_m1times(rNodes)/length(all_module_bird));

end
std_tag_of_diffRoot_t = nanstd(mean_tag_t);


h = plot(g,'MarkerSize',5,'linewidth',0.5,'NodeLabel','','NodeColor',color_m1times, ...
    'EdgeColor','r','ArrowSize',6,'EdgeAlpha',1);
layout(h,'layered','Sources',root_node)
highlight(h,root_node,'MarkerSize',10,'NodeColor','r')
text(h.XData(root_node),h.YData(root_node),num2str(root_node'),'HorizontalAlignment','center','Color','w','FontSize',8,'FontWeight','bold')

for i = 1 : length(root_node)
    highlight(h,TR{i},'EdgeColor',red_color(i,:))
end

text(mean(h.XData(root_node)),mean(h.YData(root_node))+1,'max-MS net','FontSize',12,'FontWeight','bold')
set(gca,'XTick',[],'YTick',[],'ZTick',[])
XX = gca().XLim;
xlim([XX(1)*0.0 XX(2)*0.9])
axis off


subplot('Position',sub_pos(1,:))
xyz = tracks_filt(Frame_matrix(:,frame_t),2:4);
vxyz= tracks_filt(Frame_matrix(:,frame_t),6:8);
for step=1:frame_t
    all_xyz(step,:,:) = tracks_filt(Frame_matrix(:,step),2:4);
end
vxyz = normalized_vector(vxyz);

for i=1:numBird
    plot(squeeze(all_xyz(1:frame_t,i,1)),squeeze(all_xyz(1:frame_t,i,2)),'-','color',hex2rgb('CCCCCC'),'linewidth',0.25);
    hold on;
end
box on
view(-90,90)

h = plot(g,'XData',xyz(:,1),'YData',xyz(:,2),'Marker','none','LineWidth',0.5,'EdgeColor','r',...
    'NodeLabel','','NodeColor','none','ArrowSize',6,'EdgeAlpha',1);
for i = 1 : length(root_node)
    highlight(h,TR{i},'EdgeColor',red_color(i,:))
end

scatter(xyz(:,1),xyz(:,2),50,tag_m1times,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'));
scatter(xyz(root_node,1),xyz(root_node,2),200,tag_m1times(root_node)','o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'));
%quiver3(xyz(:,1),xyz(:,2),xyz(:,3),vxyz(:,1),vxyz(:,2),vxyz(:,3),'color',[0.2 0.2 0.2],'linewidth',1,'AutoScaleFactor',0.3)
colormap(roma);
caxis([0.5 max(tag_m1times)+.5])
%hh = text(xyz(:,1),xyz(:,2),num2str([1:numBird]'),'color','w','HorizontalAlignment','center','FontSize',6);
%set(hh(root_node),'Color','w','FontWeight','bold','FontSize',10);
text(xyz(root_node,1),xyz(root_node,2),num2str(root_node'),'Color','w','FontWeight','bold','FontSize',8,'HorizontalAlignment','center');

title(['Frame = ' num2str(frame_t)])
hold off;
% axis equal
% axis(BD);
axis off


end

function [all_MS_net,aveKij_Tree,aveKij_Tree_t,S_pass2Node,corr_S_T,ave_Cr,corr_Cr_T,std_tag_of_diffRoot,std_tag_of_diffRoot_t,std_spatial_of_diffRoot,std_spatial_of_diffRoot_t] = for_a_track(forward_or_omni,consider_ego_motion,is_retina,next_tau_slice,all_MS_type,Frame_matrix,tracks_filt,turn_interval,all_module_bird_1times,all_module_bird_2times)


for i = 1 : length(all_MS_type)
    [all_MS_net{1,i},spatial_value_ij{1,i}] = Cal_1_MSnet_for_a_track(forward_or_omni,consider_ego_motion,is_retina,all_MS_type{i},Frame_matrix,tracks_filt);

    [aveKij_Tree(i),aveKij_Tree_t(i,:)] = Cal_2_Kij_for_a_track(next_tau_slice,all_MS_net{1,i},Frame_matrix,tracks_filt);
    [S_pass2Node(i,:),corr_S_T(i),ave_Cr(i,:),corr_Cr_T(i)] = Cal_3_RootNode_to_allNode_v2(all_MS_net{1,i},turn_interval);
    [std_tag_of_diffRoot(i),std_tag_of_diffRoot_t(i,:),std_spatial_of_diffRoot(i),std_spatial_of_diffRoot_t(i,:)] = Cal_45_Modules_for_a_track_v3(all_MS_net{1,i},spatial_value_ij{1,i},all_module_bird_1times,all_module_bird_2times);

end
end

