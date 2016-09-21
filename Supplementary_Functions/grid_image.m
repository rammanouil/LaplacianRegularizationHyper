
function image_a = grid_image

R = 5;
N = 75;

sampla_a = zeros(R,R^2);
sampla_a(:,1:R) = [1,0,0,0,0;
                   0,1,0,0,0;
                   0,0,1,0,0;
                   0,0,0,0,1;
                   0,0,0,1,0];
               
sampla_a(:,R+1:2*R) = 1/2*[1,0,0,0,1;
                           1,1,0,0,0;
                           0,1,1,0,0;
                           0,0,0,1,1;
                           0,0,1,1,0];
                       
sampla_a(:,2*R+1:3*R)=1/3*[1,0,0,1,1;
                           1,1,0,0,1;
                           1,1,1,0,0;
                           0,0,1,1,1;
                           0,1,1,1,0];
                       
sampla_a(:,3*R+1:4*R)=1/4*[1,0,1,1,1;
                           1,1,0,1,1;
                           1,1,1,0,1;
                           0,1,1,1,1;
                           1,1,1,1,0]; 
                       
sampla_a(:,4*R+1:5*R)=1/5*ones(R);

background_a = [0.1150,0.0741,0.2003,0.2055,0.4051]'; % as in article of TV

image_a = ones(N,N,R);
for i = 1 : N
    for j = 1 : N
        image_a(i,j,:)= background_a;
    end
end

for i = 1 : R
    for j = 1 : R
        idx_l = 5+(i-1)*15+1:5+(i-1)*15+5;
        idx_c = 5+(j-1)*15+1:5+(j-1)*15+5;
        for k = 1 : R
           for l = 1 : R
             image_a(idx_l(k),idx_c(l),:)= sampla_a(:,(i-1)*R+j);
           end
        end
    end
end

% figure, 
% for i = 1 : R
%     subplot(1,R,i);
%     imshow(image_a(:,:,i));colorbar
% end
% for i = 1 : R
%     subplot(2,round(R/2),i);
%     imshow(image_a(:,:,i)); colorbar
% end
% colormap jet
