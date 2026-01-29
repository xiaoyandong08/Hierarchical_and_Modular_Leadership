function [S_pass2Node,corr_S_T,ave_Cr,corr_Cr_T,corr_S_T1] = Cal_3_RootNode_to_allNode_v2(all_MS_net,turn_interval)


N = size(all_MS_net,2);

for t = 1 : size(all_MS_net,1)
    MS_net = squeeze(all_MS_net(t,:,:));
    [~,Cr_t(:,t)] = global_reaching_centrality(sign(MS_net));
end

for idx = 1 : N
    otherNode = [];
    for t = 1 : size(all_MS_net,1)
        MS_net = squeeze(all_MS_net(t,:,:));

        [~,Cr_t(:,t)] = global_reaching_centrality(sign(MS_net));

        TR = shortestpathtree(digraph(MS_net),idx);
        rNodes = TR.Edges.EndNodes(:,2)';

        otherNode = unique([otherNode rNodes]);
        num_pass2Node1(idx,t) = length(otherNode);
        
        num_pass2Node(idx,t) = length(rNodes);
    end
end

ave_Cr = mean(Cr_t,2)';
corr_Cr_T = corr(ave_Cr',turn_interval','type','Spearman');

S_pass2Node = sum(num_pass2Node,2)'/(size(all_MS_net,1)*(N-1));
corr_S_T = corr(S_pass2Node',turn_interval','type','Spearman');

S_pass2Node1 = sum(num_pass2Node1,2)'/(size(all_MS_net,1)*(N-1));
corr_S_T1 = corr(S_pass2Node1',turn_interval','type','Spearman');

% figure;p = plot(digraph(MS_net),'Layout','layered');highlight(p,TR,'EdgeColor','r');title(['t = ' num2str(t)])
end