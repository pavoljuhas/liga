function [d, r] = msulogo(n)
% [D, R] = MSULOGO(N), where N = 1 for "M", 2 "S", 3 "U"

if nargin < 1
    n = 1;
end

e = inline('ones(1,x)', 'x');
z = inline('zeros(1,x)', 'x');
rMSU{1} = [
    e(6), [2 3 4], 5*e(6),    z(7), 1:3,    4:5, 6*e(7);
    0:5,  [4 3 4], 5:-1:0,    0:6,  6:-1:4, 5:6, 6:-1:0;
]';

r = cat(1, rMSU{n});
d = distmx(r, 's');
if d(1) == 0
    error('repeated entry in r')
end
