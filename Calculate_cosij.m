function [all_cos_value,all_std_dist,cos_value_time] = Calculate_cosij(cos_type,Frame_matrix,tracks_filt)

plot_ij = nchoosek([1:size(Frame_matrix,1)],2);
all_cos_value = nan(size(Frame_matrix,1),size(Frame_matrix,1));
all_std_dist = nan(size(Frame_matrix,1),size(Frame_matrix,1));

cos_value_time = cell(size(Frame_matrix,1),size(Frame_matrix,1));
for k = 1 : size(plot_ij,1)
    plot_i = plot_ij(k,1);
    plot_j = plot_ij(k,2);
    cos_value = nan(1,size(Frame_matrix,2));
    dxyz = nan(1,size(Frame_matrix,2));
    for i = 1 : size(Frame_matrix,2)
        Id = Frame_matrix([plot_i plot_j],i);
        v_xyz = tracks_filt(Id,6:8);
        v_xyz = v_xyz./vecnorm(v_xyz')';
        cos_value(i) = dot(v_xyz(1,:),v_xyz(2,:));
        %cos_value(i) = norm(sum(v_xyz,1))/size(v_xyz,1);

        xyz = tracks_filt(Id,2:4);
        dxyz(i) = norm(xyz(2,:)-xyz(1,:));
    end
    cos_value_time{plot_i,plot_j} = cos_value;
    cos_value_time{plot_j,plot_i} = cos_value_time{plot_i,plot_j};
    
    if strcmp(cos_type,'mean')
        all_cos_value(plot_i,plot_j) = mean(cos_value_time{plot_i,plot_j});
    elseif strcmp(cos_type,'slope')
        temp = cos_value_time{plot_i,plot_j}(2:end)-cos_value_time{plot_i,plot_j}(1:end-1);
        all_cos_value(plot_i,plot_j) = (nanmean(temp));
    end
    
    all_cos_value(plot_j,plot_i) = all_cos_value(plot_i,plot_j);
    
    all_std_dist(plot_i,plot_j) = std(dxyz);
    all_std_dist(plot_j,plot_i) = all_std_dist(plot_i,plot_j);
end
end