function rD = mkpyramid(rABC, efg)
% MKPYRAMID     construct tetrahedron vertex using distances from base points
% rD = MKPYRAMID(rABC, efg)  where
%    rABC = [rA; rB; rC]   matrix of tetrahedron base points
%    efg = [AD, BD, CD]    is a vector of vertex distances

% line threshold
lthr = 1e-5;

% ni is unit vector in AB direction
vAB = rABC(2,:) - rABC(1,:);
c = norm(vAB);
ni = vAB/c;
% nj is a unit vector in ABC plane perpendicular to i
vAC = rABC(3,:) - rABC(1,:);
nj = vAC - dot(ni, vAC)*ni;
nj = nj/norm(nj);
% nk is a unit vector perpendicular to ABC plane
nk = cross(ni, nj);
nijk = [ni; nj; nk];

e = efg(1); f = efg(2); g = efg(3);
xA = -1/(2*c)*(c^2 + e^2 - f^2);
A = [ xA, 0, 0 ];
C = A + [dot(ni, vAC), dot(nj, vAC), 0];
xC = C(1); yC = C(2);
% vT is a translation vector from new to original coordinates
vT = rABC(1,:) - A*nijk;

h2 = e^2 - xA^2;
if sqrt(abs(h2)) < lthr
    if abs(norm(C)-g) > lthr
	return
    end
    rD = [0,0,0]*nijk + vT;
    return;
elseif h2 < 0
    rD = [];
    return;
end

yD = 1/(2*yC)*(h2 + xC^2 + yC^2 - g^2);
zD2 = h2 - yD^2;
if sqrt(abs(zD2)) < lthr
    rD = [0,yD,0]*nijk + vT;
    return;
elseif zD2 < 0
    rD = [];
    return;
end
zD = sqrt(zD2);
rD = [  [0, yD, zD]*nijk+vT;    [0, yD, -zD]*nijk+vT  ];
