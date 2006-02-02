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
% all values in C must occur somewhere in A and B
vA = sort(A(:)).'; vA(diff(vA)<=delta) = [];
vB = sort(B(:)).'; vB(diff(vB)<=delta) = [];
[ignore, idx] = sort(cij(:,1));
idx = idx(end:-1:1);
cij = cij(idx,:);
n = 1;
while n <= size(cij,1)
    CA = A(cij(1:n,2),cij(1:n,2));
    CB = B(cij(1:n,3),cij(1:n,3));
    jA = fnearbin(vA, CB(end,:));
    jB = fnearbin(vB, CA(end,:));
    if any(abs(vA(jA)-CB(end,:)) > delta) | any(abs(vB(jB)-CA(end,:)) > delta)
        cij(n,:) = [];
        continue;
    end
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
vCA = sort(CA(:)).'; vCA(diff(vCA)<=delta) = [];
vCB = sort(CB(:)).'; vCB(diff(vCB)<=delta) = [];

while 0
    nC = size(CA,1);
    xor_cnt = zeros(nC, 1);
    for i = 1:nC
        jA = fnearbin(vCA, CB(i,:));
        jB = fnearbin(vCB, CA(i,:));
        xor_cnt(i) = sum(abs(vCA(jA)-CB(i,:))>delta) + ...
            sum(abs(vCB(jB)-CA(i,:))>delta);
    end
    [m,i] = max(xor_cnt);
    if m == 0
        break
    else
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
