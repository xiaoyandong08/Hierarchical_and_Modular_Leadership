function Video_of_SupplementaryVideo2_XYZ
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

%%%%%%%%%%%


[retina_dist_ij,retina_angle_ij,spatial_value_ij,distance_ij] = Calculate_Retina_dist_of_2frame(anis_factor,Frame_matrix,tracks_filt);


forward_or_omni = 'forward';
consider_ego_motion = 'yes';
is_retina = 0;
[all_MS_net_max,ave_Gr_maxMS,Gr_t_maxMS] = Cal_1_MSnet_for_a_track(forward_or_omni,consider_ego_motion,is_retina,'max',Frame_matrix,tracks_filt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
red_color = [hex2rgb('FFCCCC');hex2rgb('F08080');hex2rgb('FF0000');hex2rgb('8B0000')];
red_color = turbo(4);
red_color(2,:) = [255 0 0]/255;
red_color(3,:) = hex2rgb('ffa500');
red_color(4,:) = hex2rgb('0072BD');
special_red_color = red_color;

% 生成Video
BD=ceil([min(tracks(:,2))-1 max(tracks(:,2))+1 ...
    min(tracks(:,3))-1 max(tracks(:,3))+1 ...
    min(tracks(:,4))-1 max(tracks(:,4))+1]);
numBird = size(Frame_matrix,1);

fun = @(m)sRGB_to_OSAUCS(m,true,true);
all_red_color = maxdistcolor(7,fun);

fig = figure;
set(gcf,'Position',[ 90    100   800   550])
for t = 1 : size(all_MS_net_max,1)

    % subplot('Position',[0.06 0.75 0.9 0.2])
    % yyaxis left; plot(Gr_t(1:t));
    % yyaxis right;plot(TT_of_Tree(1:t));
    % xlim([0 300])


    %%%%%%%%%%%%%%%%% 
    subplot('Position',[0.58 0.456 0.4 0.45])
    root_node = [];
    TR = [];
    MS_net = squeeze(all_MS_net_max(t,:,:));
    spatial_ij = squeeze(spatial_value_ij(t,:,:));
    spatial_one = nanmean(spatial_ij,1);

    [Gr_t,Cr_t,dis] = global_reaching_centrality(sign(MS_net));
    g = digraph(MS_net);
    [bins,binsize]= conncomp(g,'Type','weak');
    for i = 1 : length(binsize)
        index = find(bins==i);
        [~,index1] = max(Cr_t(index));
        root_node(i) = index(index1);
                
        temp_net = zeros(70,70);
        temp_net(index,index) = MS_net(index,index);
        new_TR{i} = digraph(temp_net);

        TR{i} = shortestpathtree(g,root_node(i));
        rNodes = TR{i}.Edges.EndNodes(:,2)';
        mean_tag_t{t,1}(i) = mean(tag_m1times(rNodes)/length(all_module_bird));
        mean_spatial_t{t,1}(i) = mean(spatial_one(rNodes));

        temp1 = rNodes(find(MS_net(rNodes,root_node(i))))';
        TR{i} = addedge(TR{i},temp1,repmat(root_node(i),length(temp1),1),ones(length(temp1),1));
    end
    std_tag_of_diffRoot_t(t,1) = std(mean_tag_t{t,1});
    std_spatial_of_diffRoot_t(t,1) = nanstd(mean_spatial_t{t,1});


    h = plot(g,'MarkerSize',5,'linewidth',1,'NodeLabel','','NodeColor',color_m1times,'EdgeColor',hex2rgb('BBBBBB'),'ArrowSize',4,'EdgeAlpha',1);
    layout(h,'layered','Sources',root_node)
    highlight(h,root_node,'MarkerSize',10,'NodeColor','r')
    text(h.XData(root_node),h.YData(root_node),num2str(root_node'),'HorizontalAlignment','center','Color','w','FontSize',8,'FontWeight','bold')
    if t==30
        red_color = special_red_color;
    else
        red_color = all_red_color(1:length(root_node),:);
    end
    for i = 1 : length(root_node)
        highlight(h,new_TR{i},'EdgeColor',red_color(i,:))
    end
    title(['max-MS network (frame = ' num2str(t) ')'],'FontSize',15.4)
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    XX = gca().XLim;
    xlim([XX(1)*0.0 XX(2)*0.9])
    axis off


    subplot('Position',[0.65 0.30 0.32 0.15])
    plot(std_tag_of_diffRoot_t,'LineWidth',2)
    xlim([0 size(all_MS_net_max,1)])
    ylabel('$\sigma^{\rm{M}}_{\rm{ST}}(t)$','Interpreter','latex','FontSize',15.4)
    set(gca,'fontsize',14)
    xlim([0 300])

    subplot('Position',[0.65 0.09 0.32 0.15])
    plot(std_spatial_of_diffRoot_t,'LineWidth',2)
    xlim([0 size(all_MS_net_max,1)])
    ylabel('$\sigma^{\rm{S}}_{\rm{ST}}(t)$','Interpreter','latex','FontSize',15.4)
    set(gca,'fontsize',14)
    xlim([0 300])
    xlabel('Frame')

    %%%%%%%%%%%%%%%%%
    subplot('Position',[0.06 0.05 0.5 0.90])
    xyz = tracks_filt(Frame_matrix(:,t),2:4);
    vxyz= tracks_filt(Frame_matrix(:,t),6:8);
    all_xyz(t,:,:) = xyz;
    vxyz = normalized_vector(vxyz);

    % for i=1:numBird
    %     plot(squeeze(all_xyz(1:t,i,1)),squeeze(all_xyz(1:t,i,2)),'-','color',hex2rgb('CCCCCC'),'linewidth',0.25);
    %     hold on;
    % end
    box on

    h = plot(g,'XData',xyz(:,1),'YData',xyz(:,2),'Marker','none','LineWidth',1,'EdgeColor','r',...
        'NodeLabel','','NodeColor','none','ArrowSize',6,'EdgeAlpha',1);
    hold on;
    for i = 1 : length(root_node)
        highlight(h,new_TR{i},'EdgeColor',red_color(i,:))
    end
    scatter(xyz(:,1),xyz(:,2),100,tag_m1times,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'));
    scatter(xyz(root_node,1),xyz(root_node,2),300,'ro','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'));
    %quiver3(xyz(:,1),xyz(:,2),xyz(:,3),vxyz(:,1),vxyz(:,2),vxyz(:,3),'color',[0.2 0.2 0.2],'linewidth',1,'AutoScaleFactor',0.3)
    colormap(roma);
    caxis([0.5 max(tag_m1times)+.5])
    hh = text(xyz(:,1),xyz(:,2),num2str([1:numBird]'),'color','w','HorizontalAlignment','center','FontSize',6);
    set(hh(root_node),'Color','w','FontWeight','bold','FontSize',10);
    
    title(['Trajectory (frame = ' num2str(t) ')'],'FontSize',15.4)
    hold off;
    axis equal
    axis(BD);
    view(-90,90)
    axis off
    %%%%%%%%%%%%%%%%%

    im(t) = getframe(fig);
end
a=VideoWriter(['Supplementary_Video_2_maxMStree_by_XYZ.mp4'],'MPEG-4');
a.FrameRate = 10;
a.Quality = 100;
open(a);
writeVideo(a,im);
close(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
