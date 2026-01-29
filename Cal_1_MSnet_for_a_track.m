function [all_MS_net,spatial_value_ij,ave_Gr,Gr_t] = Cal_1_MSnet_for_a_track(forward_or_omni,consider_ego_motion,is_retina,MS_type,Frame_matrix,tracks_filt)

anis_factor = 0;

if strcmp(consider_ego_motion,'yes')
    [retina_dist_ij,retina_angle_ij,spatial_value_ij,distance_ij] = Calculate_Retina_dist_of_2frame(anis_factor,Frame_matrix,tracks_filt);
else
    [retina_dist_ij,retina_angle_ij,spatial_value_ij,distance_ij] = Calculate_Retina_dist_of_2frame_Noego(anis_factor,Frame_matrix,tracks_filt);
end

if is_retina == 1
    MS_ij = retina_dist_ij;
else
    MS_ij = retina_angle_ij;
end

N = size(Frame_matrix,1);

for t = 1 : size(MS_ij,1)
    if strcmp(forward_or_omni,'forward')
        matrix = squeeze(MS_ij(t,:,:)).*squeeze(spatial_value_ij(t,:,:));
    elseif strcmp(forward_or_omni,'omni')
        matrix = squeeze(MS_ij(t,:,:));
    end

    if strcmp(MS_type,'max')
        [~,MS_index] = max(matrix,[],2);
    elseif strcmp(MS_type,'min')
        [~,MS_index] = min(matrix,[],2);
    elseif strcmp(MS_type,'middle')
        [~,MS_index] = sort(matrix,2,'descend');
        MS_index = MS_index(:,ceil(N/2));
    elseif strcmp(MS_type,'random')
        MS_index = randi(N,N,1);
    end

    MS_net = zeros(N,N);
    for i=1: N
        MS_net(MS_index(i),i) = 1;
    end
    all_MS_net(t,:,:) = MS_net;

    [Gr_t(t),Cr_t(t,:)] = global_reaching_centrality(sign(MS_net));

    

end

ave_Gr = mean(Gr_t);

end

