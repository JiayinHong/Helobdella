function [smoothx,smoothy] = smoothCoord(query,origin,WHATSBAD,timecourse)
%this function takes in the query coordinates and a reference origin
%coordinates as structs, and a threshold of bad jump, apply gaussian
%interpolation and return smoothed coordinates
%   query, and origin are supposed to be structs
%   try using 10 as WHATSBAD value if you have no idea what it is

% head move too much, tail seldom moves, use tail as reference
disx = query.x-origin.x;    % relative position to moving origin of tail, x coord
disy = query.y-origin.y;    % same as above, y coord
jump = sqrt([0; diff(disx).^2 + diff(disy).^2]);
figure(6);clf
plot(timecourse,disx,'k',timecourse,disy,'b')
hold on
% WHATSBAD = 10;   % hard-code, adjust according to plot of tailjump
badjump = find(jump>WHATSBAD);
plot(timecourse(badjump), 0*badjump+10, 'r*', 'MarkerSize',8);

K = length(badjump) - 1;
THR_S = 10;
for k=1:K
    if timecourse(badjump(k+1)) - timecourse(badjump(k)) < THR_S
        disx(badjump(k):badjump(k+1)-1) = nan;
        disy(badjump(k):badjump(k+1)-1) = nan;
    end
end

figure(7);clf
plot(timecourse, disx,'k', timecourse, disy, 'b')

tau = 2;
use = find(~isnan(disx));
gdisx = gaussianinterp(timecourse, timecourse(use), disx(use), tau);
gdisy = gaussianinterp(timecourse, timecourse(use), disy(use), tau);

smoothx = origin.x + gdisx;
smoothy = origin.y + gdisy;

end

