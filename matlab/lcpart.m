function [C,ij] = lcpart(A,B,delta)
% [C,IJ] = lcpart(A,B)    find largest common part of symmetric matrices A, B

% process arguments
if ~isequal(A,A.')
    error('matrix A must be symmetric');
elseif ~isequal(B,B.')
    error('matrix B must be symmetric');
end
if nargin < 3
    ndiv = double(intmax)/1000;
    delta = ( max([A(:); B(:)]) - min([A(:); B(:)]) ) / ndiv;
end

nA = size(A,1);
nB = size(B,1);

% count common numbers in every row
sA = sort(A, 2);
sB = sort(B, 2);
cij = zeros(nA*nB,3);
for i = 1:nA
    for j = 1:nB
        cij((i-1)*nB+j,:) = [sizeOfIntersect(sA(i,:), sB(j,:), delta), i, j];
    end
end
[ignore, idx] = sort(cij(:,1));
idx = idx(end:-1:1);
cij = cij(idx,:);
n = 1;
while n <= size(cij,1) & n <= cij(n,1)
    remove = find( cij(n+1:end,2) == cij(n,2) | ...
        cij(n+1:end,3) == cij(n,3)) + n;
    cij(remove,:) = [];
    n = n+1;
end
cij(n:end,:) = [];

ij = cij(:,2:3);
% now check common parts of A and B
CA = A(ij(:,1),ij(:,1));
CB = B(ij(:,2),ij(:,2));
% all values in C must occur somewhere in CA and CB
vA = sort(CA(:)).'; vA(diff(vA)<=delta) = [];
vB = sort(CB(:)).'; vB(diff(vB)<=delta) = [];

% improve with loop, check for grid27 and grid64
for i = size(CA,1):-1:1
    jA = fnearbin(vA, CB(i,:));
    jB = fnearbin(vB, CA(i,:));
    if any(abs(vA(jA)-CB(i,:))>delta) | any(abs(vB(jB)-CA(i,:))>delta)
        ij(i,:) = [];
        CA(i,:) = []; CA(:,i) = [];
        CB(i,:) = []; CB(:,i) = [];
    end
end
C = CA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function count = sizeOfIntersect(x, y, delta)
% number of equal elements in sorted vectors x, y
nx = length(x);
ny = length(y);
count = 0; j = 1;
for i = 1:nx
    while j <= ny  &  y(j) < x(i)-delta
        j = j + 1;
    end
    if j > ny;
        break
    elseif abs(y(j) - x(i)) <= delta
        count = count + 1;
        j = j+1;
    end
end
