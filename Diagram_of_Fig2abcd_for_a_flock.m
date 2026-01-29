function Diagram_of_Fig2abcd_for_a_flock
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

[all_module_bird,all_module_net] = Generate_Modularity(One_LFnet);


for i = 1 : size(all_module_bird,1)
    if length(all_module_bird{i})>5
        [module_bird_2times{i},module_net_2times{i}] = Generate_Modularity(One_LFnet(all_module_bird{i},all_module_bird{i}));
    else
        module_bird_2times{i} = {[1:1:length(all_module_bird{i})]'};
        module_net_2times{i} = 0;
    end
    for j = 1 : size(module_bird_2times{i},1)
        module_bird_2times{i}{j,2} = all_module_bird{i}(module_bird_2times{i}{j,1});
        module_bird_2times{i}{j,3} = turn_interval(module_bird_2times{i}{j,2});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_figure = 0;
node_label = [];

is_highlight = 0;

% plot 1-time module
Plot_Module_LFnet_and_Traj_v2(save_figure,node_label,Frame_matrix,tracks_filt,One_LFnet,all_module_bird,all_module_net,'group05',1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot Fig.2a
[~,Cr] = global_reaching_centrality(sign(abs(all_module_net))');
[~,sort_module1_index] = sort(Cr,'descend');

load('roma.mat');
roma = roma(ceil(linspace(1,size(roma,1),size(all_module_net,1))),:);

bird_of_module1 = zeros(1,70);
for i = 1 : length(Cr)
    bird_of_module1(all_module_bird{sort_module1_index(i)}) = i;
end

snap_frame = 190;

fig1 = figure;
set(gcf,'Position',[ 9    61   858   877])
subplot('position',[0.05 0.58 0.4 0.40]);
set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])
hold on;box on
for i = 1 : snap_frame %size(Frame_matrix,2)    
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    plot(xyz(:,1),xyz(:,2),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
end
set(gca,'fontsize',14,'XTick',[],'YTick',[],'ZTick',[])
xyz=tracks_filt(Frame_matrix(:,snap_frame),2:3);
v_xyz = normalized_vector(tracks_filt(Frame_matrix(:,snap_frame),6:7));
quiver(xyz(:,1),xyz(:,2),v_xyz(:,1),v_xyz(:,2),'color',[0.2 0.2 0.2],'linewidth',1,'AutoScaleFactor',0.7)
scatter(xyz(:,1),xyz(:,2),200,bird_of_module1,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'));
colormap(roma)
% view(-67,77)
view(-90,90)

xyz=tracks_filt(Frame_matrix(:,snap_frame),2:3);
v_xyz = normalized_vector(tracks_filt(Frame_matrix(:,snap_frame),6:7));
fig1 = figure;
set(gcf,'Position',[ 9    61   370   200])
for i = 1 : 6
    subaxis(2,3,i,'SpacingVertical',0.02,'SpacingHorizontal',0.02,'MarginLeft',.02,'MarginRight',.02, 'MarginTop',.02,'MarginBottom',.02)
    hold on;box on;
    if i==1
        bird_index = 15;
        for ii = snap_frame + 1 : size(Frame_matrix,2)
            xyz_t=tracks_filt(Frame_matrix(bird_index,ii),2:4);
            plot(xyz_t(:,1),xyz_t(:,2),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
        end
        quiver(xyz(bird_index,1),xyz(bird_index,2),4*v_xyz(bird_index,1),4*v_xyz(bird_index,2),'color',[0.2 0.2 0.2],'linewidth',0.5,'AutoScale','off')
        scatter(xyz(bird_index,1),xyz(bird_index,2),60,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'),'MarkerFaceColor',roma(1,:));

        bird_index = 8;
        for ii = snap_frame + 1 : size(Frame_matrix,2)
            xyz_t=tracks_filt(Frame_matrix(bird_index,ii),2:4);
            plot(xyz_t(:,1),xyz_t(:,2),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
        end
        quiver(xyz(bird_index,1),xyz(bird_index,2),4*v_xyz(bird_index,1),4*v_xyz(bird_index,2),'color',[0.2 0.2 0.2],'linewidth',0.5,'AutoScale','off')
        scatter(xyz(bird_index,1),xyz(bird_index,2),60,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'),'MarkerFaceColor',roma(7,:));
    else
        bird_index = all_module_bird{sort_module1_index(i)};    
        for ii = snap_frame + 1 : size(Frame_matrix,2)
            xyz_t=tracks_filt(Frame_matrix(bird_index,ii),2:4);
            plot(xyz_t(:,1),xyz_t(:,2),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
        end
        quiver(xyz(bird_index,1),xyz(bird_index,2),4*v_xyz(bird_index,1),4*v_xyz(bird_index,2),'color',[0.2 0.2 0.2],'linewidth',0.5,'AutoScale','off')
        scatter(xyz(bird_index,1),xyz(bird_index,2),60,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'),'MarkerFaceColor',roma(i,:));

    end
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    xlim([-5 20]);
    ylim([-15 15]);    
    view(-90,90)
end

end