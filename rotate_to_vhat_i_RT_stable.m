function [XvRotyz,rthetaphi,RT]=rotate_to_vhat_i_RT_stable(xyz,vhat_i,i)
% vhat_i为行向量
xyz_i=xyz(i,:);  % 行向量
% 个体i的局部坐标系，时刻：t_neg_tau
v_copy=zeros(numel(vhat_i(:,1)),2);
v_copy(:,1)=atan2(vhat_i(:, 3),sqrt(vhat_i(:,1).^2 +vhat_i(:,2).^2));  % elevation
v_copy(:,2)=atan2(vhat_i(:, 2), vhat_i(:, 1));  % azimuth
Xv_focal=xyz-xyz(i,:);

%XvRot:rotate and translate each position around the focal individual index
% 先绕z转
Rz=[cos(-v_copy(2)), -sin(-v_copy(2)),0;
    sin(-v_copy(2)), cos(-v_copy(2)),0;
    0,0,1];

Ry=[cos(-v_copy(1)), 0, -sin(-v_copy(1));
    0, 1, 0;
    sin(-v_copy(1)),0, cos(-v_copy(1))];

% 绘制坐标系验证
g_hat = [0,0,-1]';  % 重力向量（单位）
frame_x_hat = vhat_i';  % 前向为x，列向量
frame_y_hat = -cross(g_hat, frame_x_hat);  % 左方为y
frame_z_hat = cross(frame_x_hat, frame_y_hat);   % 上方为z，右手系
local_axes = [frame_x_hat,frame_y_hat,frame_z_hat];
RT=Ry * Rz * local_axes;
 
% figure(1);clf;hold on;view(3);grid on;
% quiver3(0,0,0,1,0,0,'Color',[1,0,0]);
% quiver3(0,0,0,0,1,0,'Color',[1,0,0]);
% quiver3(0,0,0,0,0,1,'Color',[1,0,0]);
% 
% quiver3(0,0,0,trans_local_axes(1,1),trans_local_axes(2,1),trans_local_axes(3,1),'Color',[1,1,0]);
% quiver3(0,0,0,trans_local_axes(1,2),trans_local_axes(2,2),trans_local_axes(3,2),'Color',[1,1,0]);
% quiver3(0,0,0,trans_local_axes(1,3),trans_local_axes(2,3),trans_local_axes(3,3),'Color',[1,1,0]);
% xlabel('x');ylabel('y');zlabel('z');


XvRotyz = Ry * Rz * Xv_focal';  % OK
XvRotyz=XvRotyz';

[az,elev,r] = cart2sph(XvRotyz(:,1),XvRotyz(:,2),XvRotyz(:,3));
rthetaphi = [r elev az];


% R_yz = Ry * Rz;
% detR_yz=det(R_yz);


end