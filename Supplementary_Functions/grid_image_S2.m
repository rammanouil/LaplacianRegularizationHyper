function image_a = grid_image_S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gererate fractional abundances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pure pixels
%x1 = eye(5); % zeros(5);x1(5,:) = ones(1,5);
x1 = zeros(15,1); x1(1,:) = 1;

% mixtures with two materials
%x2 = x1 + circshift(eye(5),[1 0]); x2 = repmat(x2(:,1),1,5);
x2 = zeros(15,1); x2(2,:) = 1; x2(3,:) = 1;

% mixtures with three materials
%x3 = x2 + circshift(eye(5),[2 0]); x3 = repmat(x3(:,1),1,5);
x3 = zeros(15,1); x3(4,:) = 1; x3(5,:) = 1;x3(6,:) = 1;

% mixtures with four  materials
%x4 = x3 + circshift(eye(5),[3 0]); x4 = repmat(x4(:,1),1,5);
x4 = zeros(15,1); x4(7,:) = 1; x4(8,:) =1;x4(9,:) = 1;x4(10,:) =1;

% mixtures with four  materials
%x5 = x4 + circshift(eye(5),[4 0]); x5 = repmat(x5(:,1),1,5);
x5 = zeros(15,1); x5(11,:) = 1; x5(12,:) = 1;x5(13,:) = 1;x5(14,:) = 1;x5(15,:) = 1;


% normalize
x2 = x2/2;
x3 = x3/3;
x4 = x4/4;
x5 = x5/5;


% background (random mixture)
%x6 = dirichlet(ones(p,1),1)';
x6 = [0 0 0 0 0 0 0 0 0 0 0.1149 0.0741  0.2003 0.2055, 0.4051]';   % as in the paper

% build a matrix
xt = [x1 x2 x3 x4 x5 x6];


% build image of indices to xt
imp = zeros(3);
imp(2,2)=1;

imind = [imp*1  imp*1 imp*1 imp*1 imp*1;
    imp*2  imp*2 imp*2 imp*2 imp*2;
    imp*3  imp*3 imp*3 imp*3 imp*3;
    imp*4  imp*4 imp*4 imp*4 imp*4;
    imp*5  imp*5 imp*5 imp*5 imp*5];

imind = kron(imind,ones(5));

% set backround index
imind(imind == 0) = 6;

% generare frectional abundances for all pixels
[nl,nc] = size(imind);
np = nl*nc;     % number of pixels
for i=1:np
    X(:,i) = xt(:,imind(i));
end

image_a = reshape(X',nl,nc,15);

% %  image endmember 1
% Xim = image_a;
% figure;
% imagesc(Xim(:,:,1))
% title('Frational abundance of endmember 1')
% colormap pink;
% colorbar; 
% 
% 
% %  image endmember 2
% figure;
% imagesc(Xim(:,:,2))
% title('Frational abundance of endmember 2')
% colormap pink;
% colorbar; 
% 
% 
% %  image endmember 3
% figure;
% imagesc(Xim(:,:,3))
% title('Frational abundance of endmember 3')
% colormap pink;
% colorbar; 
% 
% 
% %  image endmember 4
% figure;
% imagesc(Xim(:,:,4))
% title('Frational abundance of endmember 4')
% colormap pink;
% colorbar; 
% 
% 
% %  image endmember 5
% figure;
% imagesc(Xim(:,:,5))
% title('Frational abundance of endmember 5')
% colormap pink;
% colorbar; 
% 
