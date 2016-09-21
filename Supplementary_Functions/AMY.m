
function IDX = AMY(A,k)

Deg = diag(sum(A,1));
Ddemi = Deg^(-0.5);
L = Ddemi*A*Ddemi; % step 2 Normalized affinity 

[X,D] = eigs(A,k); % step 3

Y = normrow(X); % step 4

IDX = kmeans(Y,k);

% IDX_image = reshape(IDX,h,w);
% 
% figure; imshow(mat2gray(IDX_image)); colormap cool; colorbar  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
