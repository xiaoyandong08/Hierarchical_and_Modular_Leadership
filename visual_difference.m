function vis_vec=visual_difference(visual_area,vis_t,vis_t_)
% vis_t: pre
% vis_t_: next
% 1st is 俯仰角[-pi/2,pi/2]
% 2rd is 偏航角[-pi,pi]
if vis_t(2)>=-visual_area/2 & vis_t(2)<=visual_area/2
    vis_vec=vis_t_-vis_t;
    if abs(vis_vec(2))>pi
        vis_vec(2)=-1*sign(vis_vec(2))*(2*pi-abs(vis_vec(2)));
    end
    % if abs(vis_vec(1))>pi/2
    %     vis_vec(1)=-1*sign(vis_vec(1))*(pi-abs(vis_vec(1)));
    % end
else
    vis_vec = zeros(1,2);
end


end