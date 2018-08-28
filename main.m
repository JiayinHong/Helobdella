%% preprocess and store coordinates into a struct
% drag the csv file to workspace and it will be stored as a table
% change the name of the table to 'coords'
if width(coords)>=13
    coords(:,[4,7,10,13])=[];  % delete the probability columns
    tableVariableNames = {'tstamp','headX','headY','neckX','neckY','midX','midY','tailX','tailY'};
    coords.Properties.VariableNames = tableVariableNames;
end

clear str
framerate_Hz = 10;
filename = ['leech-',datestr(now,'mmmdd-yyyy'),'.mat'];     % maybe hard-coded

str.headx = coords{:,'headX'};
str.heady = coords{:,'headY'};
str.neckx = coords{:,'neckX'};
str.necky = coords{:,'neckY'};
str.midx = coords{:,'midX'};
str.midy = coords{:,'midY'};
str.tailx = coords{:,'tailX'};
str.taily = coords{:,'tailY'};

N = height(coords);
str.tt = [0:N-1]' / framerate_Hz;   % total time in seconds
save(filename, 'str')

%% load file0000 and plot head/tail trace
save_dir = '../sampleFigures/0821trial1/';
if ~isdir(save_dir)
    mkdir(save_dir)
end

close all
myrgb = imread('file02.png');   % modified in paint to remove the worm from BW background

figure; imshow(myrgb); hold on
c = linspace(1,10,N);
scatter(str.headx,str.heady,[],c,'.'); h=title('head trace raw');
cb = colorbar;  cb.Ticks(2:end-1)=[];   cb.TickLabels{1}='start';   cb.TickLabels{2}='end';
% saveas(gcf,fullfile(save_dir,[h.String,'.png']))

figure; imshow(myrgb); hold on
scatter(str.tailx,str.taily,[],c,'.'); h=title('tail trace raw');
cb = colorbar;  cb.Ticks(2:end-1)=[];   cb.TickLabels{1}='start';   cb.TickLabels{2}='end';
% saveas(gcf,fullfile(save_dir,[h.String,'.png']))

% get background BW binary
mygray = rgb2gray(myrgb);
mybinary = imbinarize(mygray,'adaptive','sensitivity',0.6);
figure; imshow(mybinary)
% saveas(gcf,fullfile(save_dir,'BlackWhite binary.png'))

%% figure showing head movement
fontsize = 16;

hx = str.headx-str.headx(1);    % relative position to moving origin of head, x coord
hy = str.heady-str.heady(1);    % same as above, y coord
headjump = sqrt([0; diff(hx).^2 + diff(hy).^2]);

figure(1);clf
set(gcf,'position',[111 569 1135 561])
subplot(2,1,1)
plot(str.tt,hx,'k',str.tt,hy,'b','LineWidth',1.5)
hold on
WHATSBAD = 30;   % hard-code, adjust according to plot of headjump
badjump = find(headjump>WHATSBAD);  % | str.headtail<str.headmid | str.headtail<str.tailmid);
plot(str.tt(badjump), 0*badjump+10, 'r*', 'MarkerSize',8);
title('head movement - raw')
legend({'x','y','artifacts'},'location','best')
set(gca,'FontSize',fontsize)
set(gca,'XTick',0:120:720)
xlim([0 721])

% then remove artifacts and apply gaussian interpolation
K = length(badjump) - 1;
THR_S = 10; % time thresh for artifacts, as they usually emerge as clustered bursts
for k=1:K
    if str.tt(badjump(k+1)) - str.tt(badjump(k)) < THR_S
        hx(badjump(k):badjump(k+1)-1) = nan;
        hy(badjump(k):badjump(k+1)-1) = nan;
    end
end

% subplot(3,1,2)
% plot(str.tt, hx,'k', str.tt, hy, 'b','LineWidth',1.5)
% title('head movement - remove artifacts')   % consider to remove this figure
% legend({'x','y'},'location','best')
% set(gca,'FontSize',fontsize)
% set(gca,'XTick',0:120:720)
% xlim([0 721])

tau = 1;
use = find(~isnan(hx));
str.ghx = gaussianinterp(str.tt, str.tt(use), hx(use), tau);
str.ghy = gaussianinterp(str.tt, str.tt(use), hy(use), tau);

subplot(2,1,2)
hp = plot(str.tt, str.ghx,'k', str.tt, str.ghy, 'b','LineWidth',1.5);
title('head movement - Gaussian smoothing')
set(gca,'FontSize',fontsize)
set(gca,'XTick',0:120:720)
xlabel('time in seconds')
xlim([0 721])
hold on

% draw gray lines representing rough surface
roughHead=[];   % record frame N.O while head dot is on rough
for i=1:N
    if isnan(str.headx(i))
        % do nothing
    elseif mybinary(round(str.heady(i)),round(str.headx(i)))==0     % caution: y coords go first!
        roughHead=[roughHead;i];
    end   
end
roughHead = roughHead/framerate_Hz;
linkaxes
h=line([roughHead roughHead],get(gca,'ylim'),'color',[0.7 0.7 0.7]);
uistack(h,'bottom')
legend([hp',h(1)],{'x','y','rough'},'location','best');
saveas(gcf,fullfile(save_dir,'head movement.png'))

%% figure showing tail movement
tx = str.tailx-str.tailx(1);    % relative position to moving origin of tail, x coord
ty = str.taily-str.taily(1);    % same as above, y coord
tailjump = sqrt([0; diff(tx).^2 + diff(ty).^2]);

figure(2);clf
set(gcf,'position',[111 569 1135 561])
subplot(2,1,1) 
plot(str.tt,tx,'k',str.tt,ty,'b','LineWidth',1.5)
hold on
WHATSBAD = 5;   % hard-code, adjust according to plot of tailjump
badjump = find(tailjump>WHATSBAD);  % | str.headtail<str.headmid | str.headtail<str.tailmid);
plot(str.tt(badjump), 0*badjump-66, 'r*', 'MarkerSize',8);
title('tail movement - raw')
legend({'x','y','artifacts'},'location','best')
set(gca,'FontSize',fontsize)
set(gca,'XTick',0:120:720)
xlim([0 721])

K = length(badjump) - 1;
THR_S = 10;
for k=1:K
    if str.tt(badjump(k+1)) - str.tt(badjump(k)) < THR_S
        tx(badjump(k):badjump(k+1)-1) = nan;
        ty(badjump(k):badjump(k+1)-1) = nan;
    end
end

% subplot(3,1,2)
% plot(str.tt, tx,'k', str.tt, ty, 'b')
% title('tail jump - remove artifacts')   % consider to remove this figure
% legend({'x','y'},'location','best')
% set(gca,'FontSize',fontsize)

tau = 1;
use = find(~isnan(tx));
str.gtx = gaussianinterp(str.tt, str.tt(use), tx(use), tau);
str.gty = gaussianinterp(str.tt, str.tt(use), ty(use), tau);

subplot(2,1,2)
hp = plot(str.tt, str.gtx,'k', str.tt, str.gty, 'b','LineWidth',1.5);
title('tail movement - Gaussian smoothing')
set(gca,'FontSize',15)
set(gca,'XTick',0:120:720)
xlabel('time in secondes')
xlim([0 721])
linkaxes
hold on

% draw gray lines representing rough surface
roughTail=[];   % record frame N.O while tail dot is on rough
for i=1:N
    if isnan(str.tailx(i))
        % do nothing
    elseif mybinary(round(str.taily(i)),round(str.tailx(i)))==0    % caution: y coords go first! 
        roughTail=[roughTail;i];
    end   
end
roughTail = roughTail/framerate_Hz;

h=line([roughTail roughTail],get(gca,'ylim'),'color',[0.7 0.7 0.7]);
uistack(h,'bottom')
legend([hp',h(1)],{'x','y','rough'},'location','best');
saveas(gcf,fullfile(save_dir,'tail movement.png'))

%% calculate head/tail velocity
% for head velocity
hvel = sqrt([0;diff(str.ghx).^2+diff(str.ghy).^2]);
hvel(hvel > 5) = nan;
use = find(~isnan(hvel));
tau = 1.5;
str.ghvel = gaussianinterp(str.tt,str.tt(use),hvel(use),tau);

% for tail velocity
tvel = sqrt([0;diff(str.gtx).^2+diff(str.gty).^2]);
% tvel(tvel > 2) = nan;
use = find(~isnan(tvel));
str.gtvel = gaussianinterp(str.tt,str.tt(use),tvel(use),tau);

%% make a movie showing head and tail velocity
if 0
    figure
    set(gcf,'position',[1000 897 997 441])
    v = VideoWriter('0821trial1.avi');
    v.FrameRate=5;
    open(v)
    tic
    for i=1:N/5
        subplot(2,1,1)
        plot(str.tt(1:i*5),str.ghvel(1:i*5),'LineWidth',1.5,'color','m');
        title('Head Velocity')
        set(gca,'FontSize',fontsize)
        set(gca,'xtick',0:120:720)
        xlim([0 721])
        ylim([0 2.5])
        hold on
        subplot(2,1,2)
        plot(str.tt(1:i*5),str.gtvel(1:i*5),'LineWidth',1.5,'color','b');
        title('Tail Velocity')
        set(gca,'FontSize',fontsize)
        set(gca,'xtick',0:120:720)
        xlabel('time in seconds')
        xlim([0 721])
        ylim([0 2.5])
        hold on
        frame=getframe(gcf);
        writeVideo(v,frame);
    end
    toc
    close(v)
end

%% get smoothed head, neck, mid and tail dot position
origin.x = str.tailx(1);
origin.y = str.taily(1);
WHATSBAD = 10;
timecourse = str.tt;

query.x = str.headx;
query.y = str.heady;
[new.headx,new.heady] = smoothCoord(query,origin,WHATSBAD,timecourse);

query.x = str.neckx;
query.y = str.necky;
[new.neckx,new.necky] = smoothCoord(query,origin,WHATSBAD,timecourse);

query.x = str.midx;
query.y = str.midy;
[new.midx,new.midy] = smoothCoord(query,origin,WHATSBAD,timecourse);

query.x = str.tailx;
query.y = str.taily;
[new.tailx,new.taily] = smoothCoord(query,origin,WHATSBAD,timecourse);

%% plot segementation length and full-body length
len(:,1) = sqrt((new.headx-new.neckx).^2+(new.heady-new.necky).^2);
len(:,2) = sqrt((new.neckx-new.midx).^2+(new.necky-new.midy).^2);
len(:,3) = sqrt((new.midx-new.tailx).^2+(new.midy-new.taily).^2);
tmp1=len(:,1);
tmp1(diff(len(:,1))>5)=nan;
use=find(~isnan(tmp1));
glen(:,1) = gaussianinterp(str.tt,str.tt(use),tmp1(use),2);
glen(:,2) = len(:,2);
tmp3=len(:,3);
tmp3(tmp3>50)=nan;
len(:,3)=tmp3;
use=find(~isnan(tmp3));
glen(:,3) = gaussianinterp(str.tt,str.tt(use),tmp3(use),2);

figure; set(gcf,'position',[1000 779 1256 559])
subplot(3,1,1); plot(str.tt,glen(:,1),'LineWidth',1.5);  title('Section I length');  set(gca,'FontSize',fontsize,'xtick',0:120:720); xlim([0 721]); ylim([0 66])
subplot(3,1,2); plot(str.tt,glen(:,2),'LineWidth',1.5);  title('Section II length'); set(gca,'FontSize',fontsize,'xtick',0:120:720); xlim([0 721]); ylim([0 66])
subplot(3,1,3); plot(str.tt,glen(:,3),'LineWidth',1.5);  title('Section III length');set(gca,'FontSize',fontsize,'xtick',0:120:720); xlim([0 721]); ylim([0 66]); xlabel('time in seconds')    
linkaxes
saveas(gcf,fullfile(save_dir,'body length by section.png'))

new.fbl = sum(glen,2);   % full-body length
figure; set(gcf,'position',[1000 808 996 310])
plot(str.tt, new.fbl,'LineWidth',1.5); h=title('full body length'); set(gca,'FontSize',fontsize,'xtick',0:120:720);xlim([0 721]);ylim([0 125])
xlabel('time in seconds')
saveas(gcf,fullfile(save_dir,[h.String,'.png']))

%% the first derivative compared with head/tail velocity
figure; set(gcf,'position',[1000 924 1256 414])

subplot(2,1,1)
dfbl=[0;diff(new.fbl)]; % the first derivative of full-body length
dfbl(dfbl>1)=nan;
use=find(~isnan(dfbl));
gdfbl=gaussianinterp(str.tt,str.tt(use),dfbl(use),2);
plot(str.tt,gdfbl,'LineWidth',1.5)
set(gca,'FontSize',fontsize,'xtick',0:120:720);
xlim([0 721]);  ylim([-1.4 1.4])
lgd=legend('$$\frac{d(Body Length)}{dt}$$','location','best');
lgd.Interpreter='latex';
%%%%%%%%%%
%  annotation('textbox',[0.2,.5,.4,.3],'LineStyle','None','String'...
% ,{'$$\frac{\mathrm{\mu Ns}}{\mathrm{\mu m}}$$'},'interpreter','latex','fontsize',11','Fontname','times') 
%%%%%%%%%%

%% plot head/tail velocity
subplot(2,1,2)
hh(1)=plot(str.tt,str.ghvel,'LineWidth',1.5,'color','m');
set(gca,'FontSize',fontsize)
hold all
hh(2)=plot(str.tt,str.gtvel,'LineWidth',1.5,'color','b'); xlabel('time in seconds')
xticks(0:120:720); xlim([0 721]); ylim([-0.5 3])

% find out tail velocity peaks
[pk,lk]=findpeaks(str.gtvel(str.gtvel>0.7),str.tt(str.gtvel>0.7),'MinPeakDistance',20);
halfpk=pk/2;    % use half-maximum to get time interval
hh(3)=patch([0 0 721 721],[-0.1 -0.5 -0.5 -0.1],[0.7 0.7 0.7],'FaceAlpha',0.5);
for i=1:length(pk)
    tmp=find(str.gtvel(lk(i)*10-100:lk(i)*10+100)>halfpk(i));
    start(i)=(tmp(1)+lk(i)*10-101)/10;
    eend(i)=(tmp(end)+lk(i)*10-101)/10;
    hh(i+3)=patch([start(i) start(i) eend(i) eend(i)],[-0.1 -0.5 -0.5 -0.1],'k');
end

% add black/gray bar below representing step/exploration
line([0 721],[-0.1 -0.1],'LineWidth',1.5,'color','k')
legend(hh(1:4),{'Head Velocity','Tail Velocity','Exploration','Step'},'location','best')
saveas(gcf,fullfile(save_dir,'BLderivative.png'))

%% bar plot showing duration of locomotion phase
stepTime = eend-start;
stepTime(1)=[]; % count since the first step
exploreTime = start(2:end)-eend(1:end-1);

figure
aver=[mean(stepTime),mean(exploreTime)];
handb=bar(aver,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
handb.BarWidth=0.3;
set(gca,'FontSize',fontsize)
hold on
err=[std(stepTime),std(exploreTime)];
errorbar(aver,err,'ko','LineWidth',1.5);

tet=text(1:length(aver),[aver(1)+10,aver(2)+40],[num2str(aver','%.2f '),[char(177);char(177)],num2str(err','% .2f sec')]...
    ,'vert','bottom','horiz','center','FontSize',fontsize);     % char(177) is plus/minus sign
box off
xticklabels({'Step time','Explore time'})
ylim([0 150]);  yticklabels([])
h=title('duration of locomotion phase');
saveas(gcf,fullfile(save_dir,[h.String,'.png']))







