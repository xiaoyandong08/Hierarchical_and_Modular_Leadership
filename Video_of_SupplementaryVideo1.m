function Video_of_SupplementaryVideo1
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

all_module_bird_2times = [];
for i = 1 : size(all_module_bird,1)
    if length(all_module_bird{i})>5
        [module_bird_2times{i},module_net_2times{i}] = Generate_Modularity(One_LFnet(all_module_bird{i},all_module_bird{i}));
    else
        module_bird_2times{i} = {[1:1:length(all_module_bird{i})]'};
        module_net_2times{i} = 0;
    end
    for j = 1 : size(module_bird_2times{i},1)
        module_bird_2times{i}{j,2} = all_module_bird{i}(module_bird_2times{i}{j,1});
        all_module_bird_2times = [all_module_bird_2times;module_bird_2times{i}(j,2)];
    end
end

[~,Cr] = global_reaching_centrality(sign(abs(all_module_net))');
[~,sort_module1_index] = sort(Cr,'descend');

load('roma.mat');
roma = roma(ceil(linspace(1,size(roma,1),size(all_module_net,1))),:);

tag_m1times = zeros(1,70);
for i = 1 : length(Cr)
    tag_m1times(all_module_bird{sort_module1_index(i)}) = i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snap_frame = 190;

xlims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))];
ylims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))];
zlims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))];

fig1 = figure;
set(gcf,'Position',[ 9    61   858   600])

for step = 1 : size(Frame_matrix,2)

    subplot('position',[0.05 0.25 0.44 0.7]);
    set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
        'XMinorTick','on','YMinorTick','on','boxstyle','full');
    box on
    Id = Frame_matrix(:,step);
    xyz=tracks_filt(Id,2:4);
    plot(xyz(:,1),xyz(:,2),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
    hold on;
    for i = 1 : step %size(Frame_matrix,2)
        Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
        xyz=tracks_filt(Id,2:4);
        plot(xyz(:,1),xyz(:,2),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
    end
    set(gca,'fontsize',14,'XTick',[],'YTick',[],'ZTick',[])
    xyz=tracks_filt(Frame_matrix(:,step),2:3);
    v_xyz = normalized_vector(tracks_filt(Frame_matrix(:,step),6:7));
    quiver(xyz(:,1),xyz(:,2),v_xyz(:,1),v_xyz(:,2),'color',[0.2 0.2 0.2],'linewidth',1,'AutoScaleFactor',0.4)
    scatter(xyz(:,1),xyz(:,2),150,tag_m1times,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'));
    colormap(roma)
    xlabel('x');ylabel('y');
    % view(-67,77)
    view(-90,90)
    hold off
    xlim(xlims)
    ylim(ylims)
    zlim(zlims)
    title(['F70, frame = ' num2str(step)],'FontSize',12)
    %%%%%%%%%%%%%%%
    
    for i = 1 : 7
        bird_index = all_module_bird{sort_module1_index(i)};

        if i<=3
            subaxis(4,6,i+3,1,'SpacingVertical',0.04,'SpacingHorizontal',0.02,'MarginLeft',.02,'MarginRight',.02, 'MarginTop',.05,'MarginBottom',.05)
        elseif i>3 & i<=6
            subaxis(4,6,i,2,'SpacingVertical',0.02,'SpacingHorizontal',0.02,'MarginLeft',.02,'MarginRight',.02, 'MarginTop',.05,'MarginBottom',.05)
        elseif i>=7
            subaxis(4,6,i-3,3,'SpacingVertical',0.02,'SpacingHorizontal',0.02,'MarginLeft',.02,'MarginRight',.02, 'MarginTop',.05,'MarginBottom',.01)

        end
        box on
        xyz=tracks_filt(Frame_matrix(bird_index,step),2:3);
        v_xyz = normalized_vector(tracks_filt(Frame_matrix(bird_index,step),6:7));

        plot(xyz(:,1),xyz(:,2),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
        hold on;
        quiver(xyz(:,1),xyz(:,2),4*v_xyz(:,1),4*v_xyz(:,2),'color',[0.2 0.2 0.2],'linewidth',0.5,'AutoScale','off')
        scatter(xyz(:,1),xyz(:,2),80,'o','filled','LineWidth',0.5,'MarkerEdgeColor',hex2rgb('262626'),'MarkerFaceColor',roma(i,:));
        colormap(roma)
        set(gca,'fontsize',14,'XTick',[],'YTick',[],'ZTick',[])
        %xlabel('x');ylabel('y');
        % view(-67,77)
        view(-90,90)
        hold off
        xlim(xlims)
        ylim(ylims)
        zlim(zlims)
        title(['M' num2str(i)],'FontSize',12)
    end
    %%%%%%%%%%%%%%%

    subplot('position',[0.05 0.07 0.93 0.1]);
    sort_Cr = sort_module1_index;
    r_agent = 0.4;
    cumsum_num_m1times = [0 cumsum(cellfun(@length,all_module_bird(sort_Cr)))'];

    v_xyz = normalized_vector(tracks_filt(Frame_matrix(cell2mat(all_module_bird_2times),step),6:7));
    quiver([1:70]',ones(70,1),-v_xyz(:,2),v_xyz(:,1),'color',[0.2 0.2 0.2],'linewidth',0.5,'AutoScaleFactor',0.3)
    hold on;

    for i = 1 : length(sort_Cr)
        arrayfun(@(x,y) rectangle('Position', [x-r_agent, y-r_agent, r_agent*2, r_agent*2], 'Curvature', [1 1], 'EdgeColor', hex2rgb('262626'), 'FaceColor', roma(i,:)), ...
            [cumsum_num_m1times(i)+1:1:cumsum_num_m1times(i+1)],1*ones(1,length(all_module_bird{sort_Cr(i)})))
        xline(cumsum_num_m1times(i+1)+0.5,'Color','r','LineWidth',1);

        if size(module_bird_2times{sort_Cr(i)},1)>1
            [~,Cr1]=global_reaching_centrality(sign(abs(module_net_2times{sort_Cr(i)}))');
            [~,sort_Cr1] = sort(Cr1,'descend');

            cumsum_num_m2times = [cumsum(cellfun(@length,module_bird_2times{sort_Cr(i)}(sort_Cr1,2)))'];
            arrayfun(@(x) xline(cumsum_num_m1times(i)+x+0.5,'Color','r','LineWidth',1),cumsum_num_m2times(1:end-1))
        end
    end
    title(['2-time modules'],'FontSize',14)
    hold off
    % ylim([-0.8 4.8]);
    % xlim([0.5 70.5])
    axis equal
    axis off

    im(step) = getframe(fig1);
end

a=VideoWriter(['Supplementary_Video_1_Modules'],'MPEG-4');
a.FrameRate = 30;
a.Quality = 100;
open(a);
writeVideo(a,im);
close(a)


end