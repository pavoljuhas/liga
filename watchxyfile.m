function watchxyfile(filename, plotstyle)
% WATCHXYFILE(FILENAME, PLOTSTYLE)

if nargin < 2
    plotstyle = 'b*';
end

if exist(filename, 'file')
    xy = rhead(filename);
else
    xy = [];
end
x = []; y = [];
if ~isempty(xy)
    x = xy(:,1); y = xy(:,2);
end
hp = plot(mean(xlim), mean(ylim), plotstyle);
set(hp, 'XData', x, 'YData', y);
axis equal;
alim = 1.05 * max(abs([1; x(:); y(:)]));
axis([-alim, alim, -alim, alim]);
mhs = clock; mhs = mhs(4:6); v = 1;

ht = title(sprintf('%02.0f:%02.0f:%02.2f version %i', mhs, v));
fileinfo = dir(filename);
hpb = findobj(gcf, 'type', 'uicontrol', 'tag', 'watchxyfilePushButton');
if isempty(hpb)
    hpb = uicontrol('Style', 'Pushbutton',...
		    'Units', 'Normalized',...
		    'Position', [.85 0 .15 .05],...
		    'String', 'Stop',...
		    'tag', 'watchxyfilePushButton');

end
set(hpb, 'UserData', 0, 'Callback', 'set(gcbo, ''UserData'', 1);')

while 1
    while isequal(fileinfo, dir(filename)) & ~get(hpb, 'UserData')
	pause(.5)
    end
    if get(hpb, 'UserData')
	delete(hpb);
	break;
    end
    fileinfo = dir(filename);
    xy = rhead(filename);
    if ~isempty(xy)
	x = xy(:,1); y = xy(:,2);
    else
	x = []; y = [];
    end
    alim = 1.05 * max(abs([1; x(:); y(:)]));
    axis([-1 1 -1 1]*alim);
    set(hp, 'XData', x, 'YData', y);
    mhs = clock; mhs = mhs(4:6); v = v+1;
    set(ht, 'String', sprintf('%02.0f:%02.0f:%02.2f version %i', mhs, v));
end
