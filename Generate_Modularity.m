function [module_intersect_bird,module_net] = Generate_Modularity(One_LFnet)


rng(123)
bp = Bipartite((abs(One_LFnet)));
bp.community.Detect();
plotFormat = PlotFormat();
plotFormat.back_color = [41,38,78]/255;
plotFormat.cell_color = 'white';
plotFormat.use_labels = true;
plotFormat.font_size = 8;
% figure;
% bp.plotter.SetPlotFormat(plotFormat);
% bp.plotter.PlotModularMatrix();
% figure;
% subplot(121);PlotWebs.PLOT_MATRIX((abs(One_LFnet)), plotFormat);
% subplot(122);PlotWebs.PLOT_NESTED_MATRIX((abs(One_LFnet)), plotFormat);

row_label_in_plot = bp.community.index_rows;%row label in plot
col_label_in_plot = bp.community.index_cols;
row_module_index = bp.community.row_modules;% module index for each bird
col_module_index = bp.community.col_modules;


% modularity trajectory
for i = 1 : max(row_module_index)
    for j = 1 : max(col_module_index)
        index_row = i;
        index_col = j;
        index1 = find(row_module_index(row_label_in_plot)==index_row);
        a1 = row_label_in_plot(index1);
        index1 = find(col_module_index(col_label_in_plot)==max(col_module_index)-index_col+1);
        a2 = col_label_in_plot(index1);

        module_One_LFnet{index_row,index_col} = One_LFnet(a1,a2);
        module_row{index_row,index_col} =  a1;
        module_col{index_row,index_col} =  a2;

        module_intersect_bird{index_row,index_col} = intersect(module_row{index_row,index_col},module_col{index_row,index_col});
        module_intersect_One_LFnet{index_row,index_col} = One_LFnet(module_intersect_bird{index_row,index_col},module_intersect_bird{index_row,index_col});

    end
end


% 下买的代码没办法处理例如3x4这种的modularity
% for i = 1 : max(row_module_index)^2
%     [index_row,index_col] = ind2sub([max(row_module_index) max(row_module_index)],i);
%     index1 = find(row_module_index(row_label_in_plot)==index_row);
%     a1 = row_label_in_plot(index1);
%     index1 = find(col_module_index(col_label_in_plot)==max(row_module_index)-index_col+1);
%     a2 = col_label_in_plot(index1);
% 
%     module_One_LFnet{index_row,index_col} = One_LFnet(a1,a2);
%     module_row{index_row,index_col} =  a1;
%     module_col{index_row,index_col} =  a2;
% 
%     module_intersect_bird{index_row,index_col} = intersect(module_row{index_row,index_col},module_col{index_row,index_col});
%     module_intersect_One_LFnet{index_row,index_col} = One_LFnet(module_intersect_bird{index_row,index_col},module_intersect_bird{index_row,index_col});
% 
%     % module_intersect_ICN{index_row,index_col} = ICN_corr(module_intersect_bird{index_row,index_col},module_intersect_bird{index_row,index_col});
% 
% end

module_intersect_bird = module_intersect_bird(:);
module_intersect_bird(cellfun(@length,module_intersect_bird)<1) = [];

aa = cellfun(@length,module_intersect_bird');
[~,index] = sort(aa,'descend');
module_intersect_bird = module_intersect_bird(index);

module_net = zeros(size(module_intersect_bird,1),size(module_intersect_bird,1));
for i = 1 : size(module_intersect_bird,1)
    for j = 1 : size(module_intersect_bird,1)
        if i~=j
            net = One_LFnet(module_intersect_bird{i},module_intersect_bird{j});
            if sum(net(:)~=0)>0.5*length(net(:))
                module_net(i,j) = mean(net(:));
            end
        end
    end
end

[~,Cr]=global_reaching_centrality(sign(abs(module_net))');
[~,sort_Cr] = sort(Cr,'descend');
module_net = module_net(sort_Cr,sort_Cr);
module_intersect_bird = module_intersect_bird(sort_Cr);

% [~,sort_column] = sort(sum(module_net~=0,1),'ascend');
% [~,sort_row]    = sort(sum(module_net~=0,2)','ascend');
% nestedness = Nestedness.NODF(sign(abs(module_net)) + eye(size(module_net,1)));
% NODF = nestedness.N;
% figure;
% subplot(121);PlotWebs.PLOT_MATRIX((abs(module_net)));
% subplot(122);PlotWebs.PLOT_NESTED_MATRIX((abs(module_net)));

end


