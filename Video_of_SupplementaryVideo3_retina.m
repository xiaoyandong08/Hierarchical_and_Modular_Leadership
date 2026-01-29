function Video_of_SupplementaryVideo3_retina
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
anis_factor = 0;

[retina_dist_ij_Noego,retina_angle_ij_Noego,spatial_value_ij_Noego,distance_ij_Noego] = Calculate_Retina_dist_of_2frame_Noego(anis_factor,Frame_matrix,tracks_filt);

forward_or_omni = 'forward';
consider_ego_motion = 'no';
is_retina = 1;
[all_MS_net_max,ave_Gr_maxMS,Gr_t_maxMS] = Cal_1_MSnet_for_a_track(forward_or_omni,consider_ego_motion,is_retina,'max',Frame_matrix,tracks_filt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : size(Frame_matrix,2)  
    v = tracks(Frame_matrix(:,i),6:8);
    xyz = tracks(Frame_matrix(:,i),2:4);
    for j = 1 : size(Frame_matrix,1)
        r{j,i} = get_visual_field_3D_without_distance(v,xyz,j);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image_dir='mobbing5_cubemap_bird';
files = dir(image_dir);
xticks_angle = [pi pi/2 0 -pi/2 -pi];
yticks_angle = [pi/2 0 -pi/2];
xticklabels_angle = {'$\pi$','$\pi/2$','$0$','$-\pi/2$','$-\pi$'};
yticklabels_angle = {'$\pi/2$','$0$','$-\pi/2$'};

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


fig = figure;
set(gcf,'Position',[ 90    100   800   550])
for t = 1 : size(all_MS_net_max,1)

    %%%%%%%%%%%%%%%%% 
    root_node = [];
    TR = [];
    MS_net = squeeze(all_MS_net_max(t,:,:));
    spatial_ij = squeeze(spatial_value_ij_Noego(t,:,:));
    spatial_one = nanmean(spatial_ij,1);

    [Gr_t,Cr_t,dis] = global_reaching_centrality(sign(MS_net));
    all_Cr_t(:,t) = Cr_t;
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

    subplot('Position',[0.5 0.55 0.4 0.42])
    h = plot(g,'MarkerSize',5,'linewidth',2,'NodeLabel','','NodeColor',color_m1times,'EdgeColor','r','ArrowSize',4,'EdgeAlpha',1);
    layout(h,'layered','Sources',root_node)
    highlight(h,root_node,'MarkerSize',10,'NodeColor','r')
    text(h.XData(root_node),h.YData(root_node),num2str(root_node'),'HorizontalAlignment','center','Color','w','FontSize',8,'FontWeight','bold')

    fun = @(m)sRGB_to_OSAUCS(m,true,true);
    all_red_color = maxdistcolor(length(root_node),fun);
    red_color = all_red_color;

    for i = 1 : length(root_node)
        highlight(h,new_TR{i},'EdgeColor',red_color(i,:))
    end
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    XX = gca().XLim;    YY = gca().YLim;
    xlim([XX(1)*0.0 XX(2)*0.9])
    text(mean(XX),YY(2),['max-MS network (frame = ' num2str(t) ')'],'FontSize',15.4,'FontWeight','bold','HorizontalAlignment','center')
    axis off

    %%%%%%%%%%%%%%%%%
    subplot('position',[0.03 0.08 0.3 0.84]);
    xyz = tracks_filt(Frame_matrix(:,t),2:4);
    vxyz= tracks_filt(Frame_matrix(:,t),6:8);
    all_xyz(t,:,:) = xyz;
    vxyz = normalized_vector(vxyz);

    % for i=1:numBird
    %     plot(squeeze(all_xyz(1:t,i,1)),squeeze(all_xyz(1:t,i,2)),'-','color',hex2rgb('CCCCCC'),'linewidth',0.25);
    %     hold on;
    % end
    box on

    h = plot(g,'XData',xyz(:,1),'YData',xyz(:,2),'Marker','none','LineWidth',2,'EdgeColor','r',...
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
    axis(BD);
    view(-90,90)
    axis off
    % axis equal

    %%%%%%%%%%%%%%%%%
    for bird_index = 24%1:7

        subplot('position',[0.38 0.05 0.6 0.5]);
        %subplot('position',pos_subplot(bird_index,:));
        %%%%%%%%%%%%%%%%%%%%
        pattern_id=sprintf('focal%d',BirdIDs(bird_index));
        pattern=sprintf('step%d.',t+1);
        filteredStruct = arrayfun(@(x) ((contains(x.name,pattern))&&(contains(x.name,pattern_id))), files, 'UniformOutput', false);
        filteredStruct_arr=cell2mat(filteredStruct);
        filteredStruct_val=files(logical(filteredStruct_arr));
        cube_img_path=sprintf('%s/%s',filteredStruct_val.folder,filteredStruct_val.name);
        cube_img=imread(cube_img_path);
        cube_img_1center=zeros(size(cube_img),'uint8');
        cube_img_1center(:,round(size(cube_img,2)/8*6+1):end,:)=cube_img(:,1:round(size(cube_img,2)/8*2),:);
        cube_img_1center(:,1:round(size(cube_img,2)/8*6),:)=cube_img(:,round(size(cube_img,2)/8*2+1):end,:);
        % cube_img_1center=cube_img_1center(:,end:-1:1,:);

        x_pixel_range = [1:1:size(cube_img_1center,2)];
        y_pixel_range = [1:1:size(cube_img_1center,1)];
        x_angle_range = linspace(pi,-pi,x_pixel_range(end));
        y_angle_range = linspace(pi/2,-pi/2,y_pixel_range(end));

        imagesc(cube_img_1center);
        axis equal
        xticks = [0, size(cube_img_1center,2)]; % 定义刻度
        yticks = [0, size(cube_img_1center,1)]; % 定义刻度
        xlim([0,size(cube_img_1center,2)]);
        ylim([0,size(cube_img_1center,1)]);

        [~,x_index] = arrayfun(@(x) min(abs(x-x_angle_range)),xticks_angle,'UniformOutput',false);
        x_index = cell2mat(x_index);
        [~,y_index] = arrayfun(@(x) min(abs(x-y_angle_range)),yticks_angle,'UniformOutput',false);
        y_index = cell2mat(y_index);
        set(gca, 'XTick',x_pixel_range(x_index),'YTick',y_pixel_range(y_index)); % 设置刻度
        set(gca, 'xticklabel', xticklabels_angle,'yticklabel', yticklabels_angle,'TickLabelInterpreter','latex'); % 设置刻度标签
        % x_labelArray = [xticklabels_angle;compose('%d',x_pixel_range(x_index));];
        % y_labelArray = [yticklabels_angle;compose('%d',y_pixel_range(y_index));];
        % set(gca, 'xticklabel', strtrim(sprintf('%s\\newline%s\n', x_labelArray{:})),'yticklabel', strtrim(sprintf('%s\\newline%s\n', y_labelArray{:}))); % 设置刻度标签
        set(gca,'FontSize',14)
        title(sprintf('focal=%d, frame=%d',bird_index,t));
        hold on;

        %%%%%%%%%%%%%%%%%%%%
        vhat_i_now =  normalized_vector(tracks_filt(Frame_matrix(bird_index,t+1),6:8));
        xyz_t_now = tracks(Frame_matrix(:,t+1),2:4);
        [xij_rel_now,R_polar_now] = rotate_to_vhat_i_RT_stable(xyz_t_now,vhat_i_now,bird_index);

        neigh_idx = [1:numBird];
        [~,x_index] = arrayfun(@(x) min(abs(x-x_angle_range)),R_polar_now(neigh_idx,3),'UniformOutput',false);
        x_index = cell2mat(x_index);
        [~,y_index] = arrayfun(@(x) min(abs(x-y_angle_range)),R_polar_now(neigh_idx,2),'UniformOutput',false);
        y_index = cell2mat(y_index);
        % scatter(x_index,y_index,50,'ro','filled')
        % text(x_index,y_index,num2str(neigh_idx'),'color','r','fontsize',10,'HorizontalAlignment','center','FontWeight','bold')
        
        %text(x_index(root_node),y_index(root_node),num2str(neigh_idx(root_node)'),'color','r','fontsize',10,'HorizontalAlignment','center','FontWeight','bold')
        for i = 1 : length(root_node)
            temp_node = new_TR{i}.Edges.EndNodes;
            temp_node = unique(temp_node(:));
            temp_node = setdiff(temp_node,bird_index);
            scatter(x_index(temp_node),y_index(temp_node),50,'o','MarkerEdgeColor',red_color(i,:),'MarkerFaceColor',red_color(i,:))
        end

        % temp_MS_net = MS_net;
        % % temp_MS_net(bird_index,:) = 0;temp_MS_net(:,bird_index) = 0;
        % temp_g = digraph(temp_MS_net);
        % h = plot(temp_g,'XData',x_index,'YData',y_index,'Marker','none','LineWidth',0.5,'EdgeColor','r',...
        %     'NodeLabel','','NodeColor','none','ArrowSize',6,'EdgeAlpha',1);
        % for i = 1 : length(root_node)
        %     highlight(h,new_TR{i},'EdgeColor',red_color(i,:))
        % end

        hold off

    end
    im(t) = getframe(fig);
end
a=VideoWriter(['Supplementary_Video_3_maxMStree_by_Retina'],'MPEG-4');
a.FrameRate = 20;
a.Quality = 100;
open(a);
writeVideo(a,im);
close(a)


end