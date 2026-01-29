function Plot_Module_LFnet_and_Traj_v2(save_figure,node_label,Frame_matrix,tracks_filt,One_LFnet,all_module_bird,all_module_net,file_tag,frame_index,cnt)

[~,sort_row]    = sort(sum(all_module_net~=0,2)','descend');
[~,sort_column] = sort(sum(all_module_net~=0,1),'descend');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 画整个的LF net
matrix = One_LFnet;
fig1 = figure;
set(gcf,'Position',[ 9    61   858   877])
ax = subplot('position',[0.05 0.05 0.5 0.5]);
hold on;box on;
imagesc(matrix(cell2mat(all_module_bird(sort_row)),cell2mat(all_module_bird(sort_column))));
set(gca,'Ydir','reverse')
colormap([parula;[1 1 1]])
% colormap([[1 1 1];parula;[1 1 1]])
axis equal
xlim([0.5 size(matrix,1)+0.5]);
ylim([0.5 size(matrix,1)+0.5]);
[rows, cols] = size(matrix);
for i = 0:cols
    line([i+0.5, i+0.5], [0.5, rows+0.5], 'Color', hex2rgb('CCCCCC'), 'LineWidth', 0.5);
end
for i = 0:rows
    line([0.5, cols+0.5], [i+0.5, i+0.5], 'Color', hex2rgb('CCCCCC'), 'LineWidth', 0.5);
end
modular_xtick = (cellfun(@length,all_module_bird(sort_column))/2)+0.5;
block_boundaries_x = cumsum(cellfun(@length,all_module_bird(sort_column))');
modular_xtick = [modular_xtick(1);modular_xtick(2:end)+block_boundaries_x(1:end-1)'];
for ii = 1 : length(sort_column)
    modular_xtick_label{ii} = ['M' num2str(sort_column(ii))];
end
set(gca,'XTick',modular_xtick','XTickLabel',modular_xtick_label,'FontSize',14)

modular_ytick = (cellfun(@length,all_module_bird(sort_row))/2)+0.5;
block_boundaries_y = cumsum(cellfun(@length,all_module_bird(sort_row))');
modular_ytick = [modular_ytick(1);modular_ytick(2:end)+block_boundaries_y(1:end-1)'];
for ii = 1 : length(sort_row)
    modular_ytick_label{ii} = ['M' num2str(sort_row(ii))];
end
set(gca,'YTick',modular_ytick','YTickLabel',modular_ytick_label,'FontSize',14)
ylabel('Module index','FontSize',14)
xlabel('Module index','FontSize',14)

for i = 1:length(all_module_bird)
    width = length(all_module_bird{i});
    bottom_x = modular_xtick(find(sort_column==i))-width/2-0.5;
    bottom_y = modular_ytick(find(sort_row==i)) - width/2-0.5;
 
    rectangle('Position', [bottom_x + 0.5, bottom_y + 0.5, width, width], 'EdgeColor', 'k', 'LineWidth', 2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 画整个的flock trajectory
subplot('position',[0.05 0.58 0.4 0.40]);
is_highlight = 0;
Plot_trajectory_for_One_Frame_onlyTraj(Frame_matrix,tracks_filt,is_highlight)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 画每个module的flock trajectory
rows = ceil(sqrt(length(all_module_bird)));
cols = ceil(length(all_module_bird)/ceil(sqrt(length(all_module_bird))));
% sub_pos = subplot_in_given_position([0.49 0.58 0.50 0.42],rows,cols);
sub_pos = subplot_in_given_position([0.52 0.58 0.45 0.40],rows,cols);

xlims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))];
ylims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))];
zlims = [floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))];
% fig2 = figure;
% set(gcf,'position',[102 481 1318 1276])
for i = 1 : length(all_module_bird)
    % subaxis(ceil(sqrt(length(all_module_bird))),ceil(length(all_module_bird)/ceil(sqrt(length(all_module_bird)))),i,...
    %     'SpacingVertical',0.02,'SpacingHorizontal',0.02)
    subplot('position',sub_pos(i,:))
    Plot_order_trajectory_for_all_module(fig1,node_label,i,...
            Frame_matrix(all_module_bird{i},:),tracks_filt,xlims,ylims,zlims);
    set(gca,'XTick',[],'YTick',[],'ZTick',[])
    %xlabel('x');ylabel('y');zlabel('z')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 画tau_ij的在module和整个网络的对比distribution, 画module network
rows = 2;
cols = 1;
sub_pos = subplot_in_given_position([0.63 0.03 0.35 0.53],rows,cols);
% subplot('position',sub_pos(1,:));PlotWebs.PLOT_NESTED_MATRIX(sign(abs(One_LFnet)));
% subplot('position',sub_pos(2,:));PlotWebs.PLOT_NESTED_MATRIX(sign(abs(all_module_net)));
subplot('position',[0.6300    0.3950    0.3500    0.1650]);
for i = 1 : size(all_module_bird,1)
    temp = One_LFnet(all_module_bird{i},all_module_bird{i});
    module_LF{i,1} = temp(temp~=0);
end
histogram(One_LFnet(One_LFnet~=0), 20, 'Normalization', 'probability', 'FaceColor', 'r', 'EdgeColor', 'n', 'FaceAlpha', 0.5);
hold on;
histogram(cell2mat(module_LF),     20, 'Normalization', 'probability', 'FaceColor', 'b', 'EdgeColor', 'n', 'FaceAlpha', 0.5);
legend('$\tau_{ij}^*$ in the LF network','$\tau_{ij}^*$ in 7 modules','Interpreter','latex','Location','best')
set(gca,'fontsize',14)
xlabel('$\tau_{ij}^*$','Interpreter','latex')
ylabel('Probability')

subplot('position',[0.58   -0.0100    0.4    0.3650]);
matrix = all_module_net;
matrix(matrix(:)==0) = inf;
min_vals = min(abs(matrix), [], 2);
mask = matrix == -min_vals;
B = matrix .* mask;
B(isnan(B)|isinf(B)) = 0;
single_modular_net = B;

G = digraph(abs(all_module_net)');
LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
p = plot(G,'layout','layered','NodeLabel',{},'LineWidth',LWidths,...
    'ArrowSize',14,'MarkerSize',30,'NodeColor',hex2rgb('F7941D'),'EdgeColor',hex2rgb('999999'),'EdgeAlpha',0.7);axis off
for ii = 1 : length(all_module_bird)
    Module_label{ii} = ['M' num2str(ii)];
end
text(p.XData,p.YData,Module_label,'Color',hex2rgb('262626'),'HorizontalAlignment','center','FontSize',14);
set(gca,'fontsize',14)
highlight(p,digraph(abs(single_modular_net)'),'EdgeColor','r')

colormap(ax,[parula;[1 1 1]])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if save_figure == 1
    saveas(gcf,['BAGE_module_large30/index' num2str(cnt) '_' file_tag '_Frame' num2str(frame_index)  '.jpg'])
    close(fig1)
end
end