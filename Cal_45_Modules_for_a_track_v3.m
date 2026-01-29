function [std_tag_of_diffRoot,std_tag_of_diffRoot_t,std_spatial_of_diffRoot,std_spatial_of_diffRoot_t,num_RootTree,ave_tree_depth] = Cal_45_Modules_for_a_track_v3(all_MS_net,all_spatial_ij,all_module_bird_1times,all_module_bird_2times)


N = size(all_MS_net,2);

tag_m1times = zeros(1,N);
for i = 1 : length(all_module_bird_1times)
    tag_m1times(all_module_bird_1times{i}) = i;
end


for t = 1 : size(all_MS_net,1)
    MS_net = squeeze(all_MS_net(t,:,:));
    spatial_ij = squeeze(all_spatial_ij(t,:,:));
    spatial_one = nanmean(spatial_ij,1);

    [Gr_t(t),Cr_t(t,:),dis] = global_reaching_centrality(sign(MS_net));
    clear tree_depth
    g = digraph(MS_net);
    [bins,binsize]= conncomp(g,'Type','weak');
    for i = 1 : length(binsize)
        index = find(bins==i);
        [~,index1] = max(Cr_t(t,index));
        component_Tree{t,1}(i) = index(index1);

        TR = shortestpathtree(g,component_Tree{t,1}(i));
        rNodes = TR.Edges.EndNodes(:,2)';
        mean_tag_t{t,1}(i) = mean(tag_m1times(rNodes)/length(all_module_bird_1times));
        
        mean_spatial_t{t,1}(i) = mean(spatial_one(rNodes));

        if length(rNodes)==0
            tree_depth(i) = 0;
        else
            tree_depth(i) = max(dis(component_Tree{t,1}(i),rNodes));
        end
    end
    tree_depth_t(t) = mean(tree_depth);
    num_RootTree_t(t) = length(component_Tree{t,1});
    std_tag_of_diffRoot_t(t) = nanstd(mean_tag_t{t,1});
    std_spatial_of_diffRoot_t(t) = nanstd(mean_spatial_t{t,1});
    
    % for i = 1 : N
    %     % TR = shortestpathtree(g,i);
    %     % rNodes = TR.Edges.EndNodes(:,2)';
    % 
    %     rNodes = find(MS_net(i,:));
    %     mean_tag_t{t,1}(i) = mean(tag_m1times(rNodes)/length(all_module_bird_1times));
    % end
    % std_tag_of_diffRoot_t(t) = nanstd(mean_tag_t{t,1});

    % g = digraph(MS_net);
    % [bins,binsize]= conncomp(g,'Type','weak');
    % for i = 1 : length(binsize)
    %     index = find(bins==i);
    %     [~,index1] = max(Cr_t(t,index));
    %     root_node(i) = index(index1);
    % end
    % 
    % figure;
    % h = plot(g,'MarkerSize',2,'linewidth',0.5,'NodeColor',hex2rgb('262626'),'ArrowSize',6,'EdgeAlpha',0.5);
    % layout(h,'layered','Sources',root_node)
    % highlight(h,root_node,'MarkerSize',10,'NodeColor','r')
    % text(h.XData(root_node),h.YData(root_node),num2str(root_node'),'HorizontalAlignment','center','Color','w','FontSize',8,'FontWeight','bold')
    % set(gca,'XTick',[],'YTick',[],'ZTick',[])
    % XX = gca().XLim;
    % xlim([XX(1)*0.0 XX(2)*0.9])
    % axis off
    % title('max-MS net')

end
ave_tree_depth = mean(tree_depth_t);
num_RootTree = mean(num_RootTree_t);
std_tag_of_diffRoot = mean(std_tag_of_diffRoot_t);
std_spatial_of_diffRoot = mean(std_spatial_of_diffRoot_t);

end