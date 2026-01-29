function Plot_trajectory_for_One_Frame_onlyTraj(Frame_matrix,tracks_filt,is_highlight)

Color = jet(size(Frame_matrix,2));
% figure
% set(gcf,'units','inches','position',[5 5.4306 5.0556 4.5694]);
box on

set(gca,'FontSize',12,'TickLength',[0.03, 0.01],...
    'XMinorTick','on','YMinorTick','on','boxstyle','full'); 
view([-162 78])
% plot3(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2),tracks_filt(Frame_matrix(Frame_matrix(:)>0),3),tracks_filt(Frame_matrix(Frame_matrix(:)>0),4),'.')
xlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),2)))])
ylim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),3)))])
zlim([floor(min(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4))) ceil(max(tracks_filt(Frame_matrix(Frame_matrix(:)>0),4)))])
hold on;box on
for i = 1 : size(Frame_matrix,2)
    
    Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
    xyz=tracks_filt(Id,2:4);
    v_xyz = tracks_filt(Id,6:8);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.','MarkerSize',2,'color',hex2rgb('CCCCCC'));
 
    % if i == size(Frame_matrix,2) 
    %     text(xyz(:,1),xyz(:,2),xyz(:,3),num2str([1:size(Frame_matrix,1)]'))
    % end
    % if i == size(Frame_matrix,2)
    %     scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
    % end
%     im(i) = getframe;
end
if is_highlight==1
    for i = 1 : size(Frame_matrix,2)
        Id = Frame_matrix(find(Frame_matrix(:,i)>0),i);
        xyz=tracks_filt(Id,2:4);
        v_xyz = tracks_filt(Id,6:8);
        plot3(xyz(24,1),xyz(24,2),xyz(24,3),'.','MarkerSize',10,'color',hex2rgb('F74461'));
    end
    Id = Frame_matrix(:,end);
    xyz=tracks_filt(Id,2:4);
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),65,[0.2 0.2 0.2],'o','filled');
    title(['Ego (red) = ' num2str(24)])
end

for i = 1 : size(Frame_matrix,1)
    Id = Frame_matrix(i,find(Frame_matrix(i,:)>0));
    if length(Id)>0
        xyz=tracks_filt(Id(end),2:4);
        scatter3(xyz(:,1),xyz(:,2),xyz(:,3),20,[0.2 0.2 0.2],'o','filled');
    end
end
set(gca,'fontsize',14)
view(-67,77)
if exist('im')
    a=VideoWriter('0.5_0.45','MPEG-4');
    a.FrameRate = 10;
    a.Quality = 100;
    open(a);
    writeVideo(a,im);
    close(a)
end
end