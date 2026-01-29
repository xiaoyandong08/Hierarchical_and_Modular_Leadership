function [retina_dist_ij,retina_angle_ij,spatial_value_ij,distance_ij] = Calculate_Retina_dist_of_2frame(Anisotropy,Frame_matrix,tracks_filt)
num_birds = size(Frame_matrix,1);

retina_angle_ij = zeros(size(Frame_matrix,2)-1,size(Frame_matrix,1),size(Frame_matrix,1));
retina_dist_ij = zeros(size(Frame_matrix,2)-1,size(Frame_matrix,1),size(Frame_matrix,1));
for i = 1 : size(Frame_matrix,1)
    for t = 1 : size(Frame_matrix,2)-1
        vhat_i_pre =  normalized_vector(tracks_filt(Frame_matrix(i,t),6:8));
        vhat_i_now =  normalized_vector(tracks_filt(Frame_matrix(i,t+1),6:8));
        xyz_t_now = tracks_filt(Frame_matrix(:,t+1),2:4);
        xyz_t_pre = tracks_filt(Frame_matrix(:,t),2:4);

        [xij_rel_now,R_polar_now] = rotate_to_vhat_i_RT_stable(xyz_t_now,vhat_i_now,i);
        [xij_rel_pre,R_polar_pre] = rotate_to_vhat_i_RT_stable(xyz_t_pre,vhat_i_now,i);


        %[R_polar_pre1,xij_rel_pre1] = get_visual_field_3D_without_distance(vhat_i_now,xyz_t_pre,1);
        %[R_polar_now1,xij_rel_now1] = get_visual_field_3D_without_distance(vhat_i_now,xyz_t_now,1);

        retina_dist = vecnorm(cell2mat(arrayfun(@(x) visual_difference(2*pi,R_polar_pre(x,2:3),R_polar_now(x,2:3)),1:num_birds,'UniformOutput',false)')');

        A1 = (1+dot(repmat([1,0,0],num_birds,1)',normalized_vector(xij_rel_now)'))/2;
        %A11 = (1+sum(repmat([1,0,0],num_birds,1).*normalized_vector(xij_rel_now),2))/2;
        retina_dist_ij(t,i,:) = retina_dist.*(A1.^Anisotropy);
        % if find(squeeze(Mij2(t,i,:))>1)
        %     find(squeeze(Mij2(t,i,:))>1)
        % end
        
        temp = dot(normalized_vector(xij_rel_now)',normalized_vector(xij_rel_pre)');
        temp(temp>1)=1;
        temp(temp<-1)=-1;
        retina_angle = acos(temp);
        % if isreal(retina_angle)~=1
        %     aa = 1;
        % end
        retina_angle_ij(t,i,:) = retina_angle.*(A1.^Anisotropy);

        spatial_value_ij(t,i,:) = A1.^1;
        distance_ij(t,i,:) = R_polar_now(:,1);
    end
end

% for t = 1 : size(retina_angle_ij,1)
%     all_retina_angle_ij{t} = squeeze(nanmean(retina_angle_ij(1:t,:,:),1));
%     all_retina_dist_ij{t}  = squeeze(nanmean(retina_dist_ij(1:t,:,:),1));
% end
end