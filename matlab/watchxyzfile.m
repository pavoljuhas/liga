function watchxyzfile(filename, max_blen)
% WATCHXYZFILE(FILENAME, MAX_BLEN)

if nargin < 2
    max_blen = 1.5;
end
atomstyle = 'b.';
bondstyle = 'k:';

org_hold = get(gca, 'NextPlot');
ha = plot3(mean(xlim), mean(ylim), mean(zlim), atomstyle);
hold on;
hb = plot3(mean(xlim), mean(ylim), mean(zlim), bondstyle);
set(gca, 'NextPlot', org_hold);

[x,y,z] = read_bga_file(filename);
set(ha, 'XData', x, 'YData', y, 'ZData', z);
[xb, yb, zb] = find_bonds(x,y,z,max_blen);
set(hb, 'XData', xb, 'YData', yb, 'ZData', zb);
axis off equal;
alim = 1.05 * max(abs([1; x; y; z]));
axis([-1 1 -1 1]*alim);
mhs = clock; mhs = mhs(4:6); v = 1;

ht = title(sprintf('%02.0f:%02.0f:%02.2f version %i', mhs, v));
fileinfo = dir(filename);
hpb = findobj(gcf, 'type', 'uicontrol', 'tag', 'watchxyzfilePushButton');
if isempty(hpb)
    hpb = uicontrol('Style', 'Pushbutton',...
		    'Units', 'Normalized',...
		    'Position', [.85 0 .15 .05],...
		    'String', 'Stop',...
		    'tag', 'watchxyzfilePushButton');

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
    [x,y,z] = read_bga_file(filename);
    set(ha, 'XData', x, 'YData', y, 'ZData', z);
    [xb, yb, zb] = find_bonds(x,y,z,max_blen);
    set(hb, 'XData', xb, 'YData', yb, 'ZData', zb);
    alim = 1.05 * max(abs([1; x; y; z]));
    axis([-1 1 -1 1]*alim);
    mhs = clock; mhs = mhs(4:6); v = v+1;
    set(ht, 'String', sprintf('%02.0f:%02.0f:%02.2f version %i', mhs, v));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y,z] = read_bga_file(filename)
% this can read XYZ and AtomEye files
x = []; y = []; z = [];
if ~exist(filename, 'file')
    return;
end
fid = fopen(filename);
[r,header] = rhead(fid, 0);
if ~isempty(strfind(header, 'Number of particles'))
    % it is atomeye file, needs to be read again
    pat   = 'H0(1,1) =';
    start = strfind(header, pat);
    a = sscanf(header(start(1)+length(pat):end), '%f', 1);
    r = rhead(fid, 0);
    r = r(:,1:3)*a - a/2.0;
end
fclose(fid);
if size(r,1)>0
    x = r(:,1); y = r(:,2); z = r(:,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xb, yb, zb] = find_bonds(x,y,z,max_blen)
dm = distmx([x(:),y(:),z(:)]);
[ib,jb] = find(dm < max_blen & dm > 0);
kb = find(ib<jb);
ib = ib(kb); jb = jb(kb);
xb = []; yb = []; zb = [];
if isempty(ib)
    return
end
xb = ones(3,length(ib))*NaN; yb = xb; zb = xb;
xb(1:3:end) = x(ib); xb(2:3:end) = x(jb);
yb(1:3:end) = y(ib); yb(2:3:end) = y(jb);
zb(1:3:end) = z(ib); zb(2:3:end) = z(jb);
xb = xb(:); yb = yb(:); zb = zb(:);
