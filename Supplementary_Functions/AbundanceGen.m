function a = AbundanceGen(R,N,rho)
% Generate abundance vectors
% R: Endmember number (Length of an abundance vector)
% N: Number of pixels
% rho: Theshold for highly mixed case. If there is a value larger than rho,
%      this abundance vector will be replaced by 1/R*ones(R,1); rho=1 by
%      defaut, meaning no abundance vector will be replaced.
% a: Output abundance matrix, with dimension R*N.

if nargin < 3
    rho=1;
end

a = rand(R,N);
a = a./ (ones(R,1)*sum(a));
idx = find(max(a)>rho);
if ~isempty(idx)
    a(:,idx) = 1/R;
end

end