
function X = normrow(M)

tmp = M.^2;
tmp = sum(tmp,2);

[n c] = size(M);

X = zeros(n,c); 

for i=1:n
    X(i,:) = M(i,:)./sqrt(tmp(i));
end