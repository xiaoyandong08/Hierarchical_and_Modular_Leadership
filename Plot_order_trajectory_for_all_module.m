function Plot_order_trajectory_for_all_module(fig,node_label,index,Frame_matrix,tracks_filt,xlims,ylims,zlims)


% xlabel('x');ylabel('y');zlabel('z')
set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
view([-162 78])
view(-67,77)

xlim(xlims)
ylim(ylims)
zlim(zlims)
% xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
% ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
% zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])

hold on;box on
for i = 1 : size(Frame_matrix,2)
    
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    %plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','color',Color(i,:));
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),2,hex2rgb('CCCCCC'),'filled');
    
    if isempty(node_label)==false & i==size(Frame_matrix,2)
        text(xyz(:,1),xyz(:,2),xyz(:,3),num2str(node_label{index}))
    end
%     im(i) = getframe;
end
for i = 1 : size(Frame_matrix,1)
    Id = Frame_matrix(i,find(Frame_matrix(i,:)>0));
    if length(Id)>0
        xyz=tracks_filt(Id(end),2:4);
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),40,[0.2 0.2 0.2],'o','filled');
    end
end
set(gca,'fontsize',14)

end