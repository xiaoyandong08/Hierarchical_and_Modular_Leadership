function Diagram_of_Fig3abcd_for_a_flock
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

[all_cos_value,~,cos_value_time] = Calculate_cosij('mean',Frame_matrix,tracks_filt);

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

load('roma.mat');
roma = roma(ceil(linspace(1,size(roma,1),size(all_module_net,1))),:);

tag_m1times = zeros(1,length(Rank));
for i = 1 : length(all_module_bird)
    tag_m1times(all_module_bird{i}) = i;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot Fig.3a
[~,Cr]=global_reaching_centrality(sign(abs(all_module_net))');
[~,sort_Cr] = sort(Cr,'descend');

r_agent = 0.4;
cumsum_num_m1times = [0 cumsum(cellfun(@length,all_module_bird(sort_Cr)))'];

figure;
set(gcf,'Position',[ 9    300   1030   100])
hold on;
for i = 1 : length(sort_Cr)
    arrayfun(@(x,y) rectangle('Position', [x-r_agent, y-r_agent, r_agent*2, r_agent*2], 'Curvature', [1 1], 'EdgeColor', hex2rgb('262626'), 'FaceColor', roma(i,:)), [cumsum_num_m1times(i)+1:1:cumsum_num_m1times(i+1)],1*ones(1,length(all_module_bird{sort_Cr(i)})))
    xline(cumsum_num_m1times(i+1)+0.5);
    text(mean([cumsum_num_m1times(i)+1 cumsum_num_m1times(i+1)]),2.7,['M' num2str(sort_Cr(i))],'HorizontalAlignment','center','FontSize',14)
end
axis equal

figure;
set(gcf,'Position',[ 9    100   1030   100])
hold on;
for i = 1 : length(sort_Cr)
    arrayfun(@(x,y) rectangle('Position', [x-r_agent, y-r_agent, r_agent*2, r_agent*2], 'Curvature', [1 1], 'EdgeColor', hex2rgb('262626'), 'FaceColor', roma(i,:)), ...
        [cumsum_num_m1times(i)+1:1:cumsum_num_m1times(i+1)],1*ones(1,length(all_module_bird{sort_Cr(i)})))
    xline(cumsum_num_m1times(i+1)+0.5);

    if size(module_bird_2times{sort_Cr(i)},1)>1
        [~,Cr1]=global_reaching_centrality(sign(abs(module_net_2times{sort_Cr(i)}))');
        [~,sort_Cr1] = sort(Cr1,'descend');

        cumsum_num_m2times = [cumsum(cellfun(@length,module_bird_2times{sort_Cr(i)}(sort_Cr1,2)))'];
        arrayfun(@(x) xline(cumsum_num_m1times(i)+x+0.5,'Color','r'),cumsum_num_m2times(1:end-1))
    end
end
axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot Fig.3b
index = 5;
figure;
hold on;box on
scatter(turn_interval,all_cos_value(index,:),30,tag_m1times,'o','filled')
colormap(roma)
scatter(turn_interval([5 24 31 48]),all_cos_value(index,[5 24 31 48]),80,'ro','filled')
xline(turn_interval(index),'-',['turning delay of bird-' num2str(index)],'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','FontSize',14)
set(gca,'fontsize',14)
xlabel('Turning delay T_{\itj} of each bird (s)')
ylabel(['$\left < {\bf{v}}_{' num2str(index) '}\cdot {\bf{v}}_{j}\right >_t$'],'Interpreter','latex')
set(gcf,'Position',[840   310   300   229])

% plot Fig.3c
index = 11;
module_index = find(cellfun(@(x) sum(find(x==index)),all_module_bird_2times)~=0);
figure;
hold on;box on
scatter(turn_interval,all_cos_value(index,:),30,tag_m1times,'o','filled')
colormap(roma)
scatter(turn_interval(all_module_bird_2times{module_index}),all_cos_value(index,all_module_bird_2times{module_index}),80,'ro','filled')
xline(turn_interval(index),'-',['turning delay of bird-' num2str(index)],'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','FontSize',14)
set(gca,'fontsize',14)
xlabel('Turning delay T_{\itj} of each bird (s)')
ylabel(['$\left < {\bf{v}}_{' num2str(index) '}\cdot {\bf{v}}_{j}\right >_t$'],'Interpreter','latex')
set(gcf,'Position',[840   310   300   229])

% plot Fig.3d

for i = 1 : size(all_cos_value,1)
    index = cellfun(@length,cellfun(@ (x) find(x==i),all_module_bird_2times,'UniformOutput',false));
    index = all_module_bird_2times{find(index)};
    aveO_inModule{i,1} = all_cos_value(i,index)';
    index = setdiff([1:size(all_cos_value,1)],index);
    aveO_outModule{i,1} = all_cos_value(i,index)'; 
end

X1 = cell2mat(aveO_inModule);
X2 = cell2mat(aveO_outModule);

figure;
bx = boxplot([X1;X2],[ones(length(X1),1);2*ones(length(X2),1);],...
    'Labels',{'within','between'});
colors = [hex2rgb('F74461');hex2rgb('1AA85A')];
h = findobj(gca,'Tag','Box');
for j = length(h):-1:1
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',0.3);
end
set(gca,'fontsize',14)
ylabel('$\left<\hat{\bf{v}}_i\cdot \hat{\bf{v}}_j\right>_t$','Interpreter','latex')
% ylim([0.9 1.0])
xlim([0.7 2.3])
set(gcf,'Position',[585   310   185   210])






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end