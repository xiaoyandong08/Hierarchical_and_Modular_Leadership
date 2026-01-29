function [aveKij_Tree,aveKij_Tree_t] = Cal_2_Kij_for_a_track(next_tau_slice,all_MS_net,Frame_matrix,tracks_filt)


N = size(all_MS_net,2);

for t = 1 : size(all_MS_net,1)
    
    MS_net = squeeze(all_MS_net(t,:,:));

    if t+next_tau_slice <= size(all_MS_net,1)
        Kij_t = Calculate_cosij('slope',Frame_matrix(:,[t+1:t+next_tau_slice]),tracks_filt);
        
        [sort_Kij,sort_index] = sort(Kij_t,2,'descend');
        
        all_agent = zeros(N,N);
        for i=1:N
            all_agent(i,sort_index(i,:)) = [1:N];
        end
        
        %present_index = sign(MS_net+MS_net');
        present_index = MS_net';
        aveKij_Tree_t(t) = nanmean(N+1-(all_agent(present_index(:)==1)-1))/N;
        aveKij_notTree_t(t) = nanmean(N+1-(all_agent(present_index(:)==0)-1))/N;

        aveKij_Tree_t1(t) = nanmean((all_agent(present_index(:)==1)-1))/N;
        aveKij_notTree_t1(t) = nanmean((all_agent(present_index(:)==0)-1))/N;

    end

end
aveKij_Tree    = nanmean(aveKij_Tree_t);

end