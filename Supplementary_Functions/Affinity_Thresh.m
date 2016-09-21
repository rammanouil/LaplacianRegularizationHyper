
function [A, xy] = Affinity_Thresh(h,N,r,dmin)

A = zeros(N,N); 
D = zeros(N,N); 

for i =1:N
    ri = r(:,i);
    D(i,i) = inf;
    
    for j = i+1:N
        rj = r(:,j);            
        D(i,j) = sum((ri - rj).^2); % +sum((rsi - rj).^2)*flag +sum((ri - rsj).^2)*flag ;
        D(j,i) = D(i,j);
    end
    
end

for i=1:N
    [tmp, idx] = sort(D(i,:)); % figure; plot(tmp)
    k = sum(tmp<dmin); %%%%%%%%% 
    A(i,idx(1:k)) = 1;
    A(idx(1:k),i) = 1;
end

% figure; imshow(mat2gray(A)); 

xy = zeros(N,2);

for i=1:N
    [l, c] = image_coordinates(h,i);
    xy(i,:) = [h+1-l, c];
end