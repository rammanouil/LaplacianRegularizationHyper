function [r, a] = hypermix(M, N, par, a)
% Hyperspectral mixing function
% M - endmember matrix
% N - number of pixels to generate
% par - mixture model (1: linear, 2: bilinear, 3, initimate)
% a - abundances (If not given, random values will be used)
% Jie CHEN. Lab. Fizeau, Unice. 2011.

[L,R] = size(M);

if nargin < 4
    a = rand(R,N);
    a = a./(repmat(sum(a(1:R,:)),R,1));
end


r = zeros(L, N);

if par == 1             % linear mixing

    for i = 1 : N
        for j = 1 : R
            r(:,i) = r(:,i) + M(:,j) * a(j,i);
        end
    end

elseif par == 2         % normal bilinear mixing

    r = hypermix(M, N, 1, a);
    for i = 1 : N
        for j = 1 : R-1
            for k = j+1 : R
                r(:, i) = r(:, i) + a(j,i)*a(k,i)* M(:,j).*M(:,k);
            end
        end
    end

elseif par == 3        % intimate mixture
    r1 = hypermix(M, N, 1, a);
    r = r1.^(0.7);
end

end