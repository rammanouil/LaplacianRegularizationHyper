
% clear ; n = 2; R = 2; s = 1; 
% [image_a] = grid_image2(n,R,s);
% figure ; imshow(image_a(:,:,1))

function Permuting_matrix = permute_grid_image(s,n,R)

h = R*n + 2*s*R ; 
N = h*h ; 
w = h;

Permuting_matrix = zeros(N,N);
background = 1:N;
background_remove = zeros(1,n*n*R*R);
count = 1;

for j=1:R % j before i to change order of squares 
    for i=1:R
        
        idx_l = s+(i-1)*(n+2*s)+1:s+(i-1)*(2*s+n)+n;
        idx_c = s+(j-1)*(n+2*s)+1:s+(j-1)*(2*s+n)+n;
        
        for k=1:n
            for l=1:n
                nn = (idx_l(k)-1)*w + idx_c(l);
                Permuting_matrix(nn,count) = 1;
                background_remove(count) = nn;
                count = count+1;
            end
        end
        
    end
end

background(background_remove) = [];

for i = 1 : length(background)
    Permuting_matrix(background(i),count) = 1;
    count = count+1;
end

% a = reshape(image_a, N, R).';
% ap = a*Permuting_matrix;
% r_l = 0.1;
% [A,~] = complete_rbf_weighted(h,N,a,r_l); 
% for i=1:N; A(i,i) = 1 ; end
% figure ; imshow(mat2gray(A));
% [Ap,~] = complete_rbf_weighted(h,N,ap,r_l); 
% for i=1:N; Ap(i,i) = 1 ; end
% figure ; imshow(mat2gray(Ap));
