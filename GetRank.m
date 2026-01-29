function [Rank,Delay,One_LFnet,DelayM,max_u]=GetRank(tracks,tau_threshold,correlation_threshold)
%%%%%%%%%%%%% function
% Rank birds who turns first, based on velocity correlation
% input: tracks (1:n,1:12) where n is number of birds
%                  (:,1) - bird id
%                  (:,2:4) - bird position
%                  (:,5) - time
%                  (:,6:8) - velocity
%                  (:,9:11) - acceleration
%                  (:,12) - wingbeat
% Associate each bird ID with a 'Rank' number
% Rank = 1:N;
% Rank = 1 is the first bird start to turn,
%      = N is the last bird turn;
% Delay = 0 for Rank = 1;
% Delay(Rank=i) = 1/(i-1) * sum ( Delay(Rank=j) + tau_ij ), where j<i
%
% Reference: ?Attanasi et al.,?Information transfer and behavioural inertia in starling
% flocks, 2014, Nature Physics
%%%%%%%%%%%%% start to calculate
Rank=[];
Delay=[];
BirdIDs=unique(tracks(:,1));
Score(1:length(BirdIDs))=0; % Socre = min: turn first
T=unique(tracks(:,5));
T=unique(tracks(:,5));
dt=T(2)-T(1); % time step
for i=1:length(BirdIDs)
    %%% for each bird, correlated with the rest of the group, delay maxtrix
    id=find(tracks(:,1)==BirdIDs(i));
    u1=tracks(id,6:8);
    xyz1 = tracks(id,2:4);
    for j=1:length(BirdIDs)
        if j==i
            DelayM(i,j)=0; % mutual delay
            continue;
        end
        id=find(tracks(:,1)==BirdIDs(j));
        u2=tracks(id,6:8);
        xyz2 = tracks(id,2:4);
        clear cor_ut cor_u
        [cor_ut,cor_u,delay_temp,max_u(i,j)]=getCor_fullLength(u1,u2,dt,tau_threshold,correlation_threshold);% get delay time by velocity correlation
        
        % figure;hold on;
        % plot(xyz1(:,1),xyz1(:,2),'r-');
        % quiver(xyz1(:,1),xyz1(:,2),u1(:,1),u1(:,2))
        % plot(xyz2(:,1),xyz2(:,2),'b-')
        % quiver(xyz2(:,1),xyz2(:,2),u2(:,1),u2(:,2))

        % [cor_ut,cor_u,delay_temp]=getCor(u1,u2,dt);
        DelayM(i,j)=delay_temp; % mutual delay
        if DelayM(i,j)<0
            Score(i)=Score(i)-1;
        elseif DelayM(i,j)>0
            Score(i)=Score(i)+1;
        end
    end
end
%%%%%%%%%%% use all birds to define the delay time -- do not use it
% Delay=mean(DelayM,2);
% Delay=Delay-min(Delay); %% delay time
% Rank=Score;

%%%%%%%%%%% only use birds turn earlier
[Score_new,I]=sort(Score);
for i=1:length(I)
    Rank(I(i))=i; % = 1 is the first bird start to turn
end
%%%% calculate the delay time
Delay(Rank==1)=0;
for i=2:length(I)
    DelaySum=0;
    idx=find(Rank==i);
    for j=1:i-1
        idy=find(Rank==j);
        Delay_temp=DelayM(idx,idy); % mutual delay between rank i and j
        DelaySum=DelaySum+Delay_temp+Delay(Rank==j);
    end
    Delay(Rank==i)=DelaySum/(i-1);
end

%%%% redo for some cases Delay < 0
id=find(Delay<0);
Delay(id)=0;
Rank(id)=1;

% %%% confrim tij=tik+tkj
% figure('units','inches','position',[2 2 3.5 3]);
% tij=[];
% tikkj=[];
% for i=1:size(DelayM,1)
%     for j=1:size(DelayM,2)
%         tij=[tij;DelayM(i,j)];
%         tikkj_temp=[];
%         for k=1:size(DelayM,1)
%             tikkj_temp=[tikkj_temp;DelayM(i,k)+DelayM(k,j)];
%         end
%         tikkj=[tikkj;mean(tikkj_temp)];
%     end
% end
% plot(tij,tikkj,'ro','MarkerFaceColor','w','MarkerSize',5);hold on
% plot([-2.5 2.5],[-2.5 2.5],'k-','LineWidth',1.5)
% axis([-2.5 2.5 -2.5 2.5]);box off
% set(gca,'Position',[0.15 0.15 0.75 0.75])
% set(gca,'FontSize',12,'FontName','times','TickLength',[0.03, 0.01],...
%     'XMinorTick','on','YMinorTick','on');
DelayM = -DelayM;

DelayM1 = DelayM;
DelayM1(DelayM1>0) = 0;
One_LFnet = DelayM1;

end

function [cor_ut,cor_u,delay,max_u]=getCor_fullLength(u1,u2,dt,tau_threshold,correlation_threshold)
%%%% get delay time based on correlation in velocity,
%%%% delay time is u1 relative to u2, i.e., delay>0, u1 turn later,
%%%% verse verse
% % %     %%%% normalize the u1 and u2
% % %     U1=sum(u1.*u1,2).^0.5;
% % %     U2=sum(u2.*u2,2).^0.5;
% % %     u1=u1./U1;
% % %     u2=u2./U2;
%%%%%%%%%%% start the calculation
%     Mid=0;
%     Slab=299;
%     u2_t=u2(Mid-Slab:Mid+Slab,:);
%     cor_ut=Slab+1-Mid:size(u1,1)-Mid-Slab-1;
all_t = [1:1:size(u1,1)];
cor_ut = [-size(u1,1)+1:1:size(u1,1)-1];
cor_u=[];
for k=1:length(cor_ut)
    index = all_t((all_t-cor_ut(k))>0&(all_t-cor_ut(k))<=size(u1,1));

    u1_t=u1(index,:);
    u2_t=u2(index-cor_ut(k),:);

    u1_t = u1_t./vecnorm((eps+u1_t)')';
    u2_t = u2_t./vecnorm((eps+u2_t)')';

    cor_u(k)=mean(sum(u1_t.*u2_t,2));

end
cor_ut=cor_ut*dt; %%% normalized to second
cor_u(cor_u<correlation_threshold)=nan;
cor_u(1:ceil(tau_threshold*length(cor_ut))) = nan;
cor_u(end-ceil(tau_threshold*length(cor_ut))+1:end) = nan;

% [I,I1] = max(cor_u);
% if length(I)>0
%     delay = cor_ut(I1);
%     max_u = I;
% else
%     delay = 0;
%     max_u = 0;
% end

[I,~] = findpeaks(cor_u);
if length(I)>0
    I = I(I == max(I));
    delay=cor_ut(cor_u==I);
    max_u = I;
else
    delay = 0;
    max_u = 0;
end
end

function [cor_ut,cor_u,delay]=getCor(u1,u2,dt)
%%%% get delay time based on correlation in velocity,
%%%% delay time is u1 relative to u2, i.e., delay>0, u1 turn later,
%%%% verse verse
% % %     %%%% normalize the u1 and u2
% % %     U1=sum(u1.*u1,2).^0.5;
% % %     U2=sum(u2.*u2,2).^0.5;
% % %     u1=u1./U1;
% % %     u2=u2./U2;
%%%%%%%%%%% start the calculation
Mid=round(size(u1,1)/2);
Slab=30;
u2_t=u2(Mid-Slab:Mid+Slab,:);
cor_ut=Slab+1-Mid:size(u1,1)-Mid-Slab-1;
cor_u=[];
for k=1:length(cor_ut)
    u1_t=u1(Mid+cor_ut(k)-Slab:Mid+cor_ut(k)+Slab,:);
    cor_u(k)=sum(sum(u1_t.*u2_t));
    cor_un=sum(sum(u1_t.^2))^0.5*sum(sum(u2_t.^2))^0.5;
    cor_u(k)=cor_u(k)/cor_un;
end
cor_ut=cor_ut*dt; %%% normalized to second
delay=cor_ut(cor_u==max(cor_u));
end

