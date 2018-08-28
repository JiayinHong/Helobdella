%% this script is used to illustrate donut-like exploration area

donutr = min(new.fbl);  % inner cycle radius
donutR = max(new.fbl);  % outter cycle radius

%% make rough-centered donut
figure
set(gcf,'position',[360 307 420 391])
xlim([0 350]); ylim([0 350]);
for x=0:140:280
    for y=0:140:280
        rectangle('Position',[x y 70 70], 'FaceColor', [0.7 0.7 0.7])
    end
end
hold all
for x=70:140:210
    for y=70:140:210
        rectangle('Position',[x y 70 70], 'FaceColor', [0.7 0.7 0.7])
    end
end
viscircles([175 175],donutr,'color','b','LineWidth',1.5,'LineStyle','-');
viscircles([175 175],donutR,'color','m','LineWidth',1.5,'LineStyle','-');
% viscircles([175 175],sqrt(105^2+30^2),'color','m','LineWidth',1.5,'LineStyle','-')

%% make smooth-centered donut
figure
set(gcf,'position',[360 307 420 391])
xlim([0 350]); ylim([0 350]);
for x=70:140:210
    for y=0:140:280
        rectangle('Position',[x y 70 70], 'FaceColor', [0.7 0.7 0.7])
    end
end
hold all
for x=0:140:280
    for y=70:140:210
        rectangle('Position',[x y 70 70], 'FaceColor', [0.7 0.7 0.7])
    end
end
viscircles([175 175],donutr,'color','b','LineWidth',1.5,'LineStyle','-');
viscircles([175 175],donutR,'color','m','LineWidth',1.5,'LineStyle','-');
