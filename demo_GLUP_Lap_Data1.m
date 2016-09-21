
% This demo reproduces the results reported in the paper "A Graph Laplacian regularization for Hyperspectral data unmixing" - ICASSP2015
% for the case of Data1 ; SNR 30 dB ;
% It solves the following optimization problem:
% 
%   min_{A}    1/2*|| S - R*A||_{F}^2 + lambda * tr(A*Lap*A^T) + mu sum_{k=1}^N || a_{k}||_2
% {subject to} A_{ij} > 0  for all  i,j  
% 			 sum_{i=1}^N A_{ij} = 1 for all j
% 
% To solve this problem we use the ADMM, and consider the following variable splitting:  
% 
%  min_{X,Y,Z} 1/2*||S-R*X||_{F}^2 + lambda*tr(Y*Lap*Y^T) + mu sum_{k=1}^N ||z_{k}||_2 + I(Z)
% {subject to} 1^*X = 1
%              X = Z 
%              X = Y
% 
% We refer to this problem as GLUP_Lap which stands for Group Lasse with Unit Sum, Positivity constraint and Laplacian regularization
% The following demo allows to compare the results of GLUP-Lap wih those obtained using FCLS and SUnSAL-TV using the data set described in 
% M.-D. Iordache, J. Bioucas-Dias, and A. Plaza, "Total variation spatial regularization for sparse hyperspectral unmixing", IEEE Transactions on 
% Geoscience and Remote Sensing, vol. PP, no. 99, pp. 1-19, 2012.
% The code for SUnSAL-TV and the dictionary were the same used for the previous paper, a demo "demo_sunsal_TV" is available at http://www.lx.it.pt/~bioucas/publications.html 
% parts of this demo are used in here (Loading dictionary - SUnSAL-TV) and
% the integrality of the demo is provided in Demo_SUnSAL_TV

clear all ; 
close all ; 
clc ; 

t_demo = clock;

%%
% Generate abundances for Data1

image_X = grid_image ; % generate Data1 used in experiments
[h, w, P] = size(image_X); % h hight w width and P nbr of endmembers 
N = w*h; % number of pixels 
X = reshape(image_X, N, P).'; % reshape abundances into P x N matrix

% figure; % show fractional abundances in aimage
% for i=1:P
%     subplot(3,5,i); imshow(image_X(:,:,i));
% end
% colormap pink; 

%%
% buid the dictionary and select endmembers

load USGS_1995_Library.mat
%  order bands by increasing wavelength
[dummy, index] = sort(datalib(:,1));
Dict =  datalib(index,4:end);

% prune the library 
% min angle (in degres) between any two signatures 
% the larger min_angle the easier is the sparse regression problem
min_angle = 4.44;       
Dict = prune_library(Dict,min_angle); % 240  signature 

% order  the columns of A by decreasing angles 
[Dict, ~, ~] = sort_library_by_angle(Dict);

% select P endmembers  from A
supp = 1:P;
M = Dict(:,supp);
[L,P] = size(M);  % L = number of bands; P = number of material

%% 
% Create data from fractional abundances X and dictionary Dict, then add noise 

SNR = 30; % [10 20 30 40]
std_noise = sqrt(sum(sum((M*X).^2))/N/L/10^(SNR/10));
noise = std_noise*randn(L,N);
Y = M*X + noise;

% create  true X wrt to the library Dict
n = size(Dict,2);
XT = zeros(n,N);
XT(supp,:) = X;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fully constrained least squares FCLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(strcat(sprintf('Waiting for fcls ...'),'\n'));

lambda = 0;
t = clock ; 
[X_hat_fcls] =  sunsal(Dict,Y,'lambda',lambda,'ADDONE','yes','POSITIVITY','yes', ...
                    'TOL',1e-4, 'AL_iters',2000,'verbose','no');
t_fcls = etime(clock,t);

SRE_fcls = 20*log10(norm(XT,'fro')/norm(X_hat_fcls-XT,'fro')); 
RMSE_fcls = Compute_RMSE(XT,X_hat_fcls);
fprintf(strcat(sprintf('FCLS : (SRE - RMSE) = (%2.3f - %2.4f)',SRE_fcls,RMSE_fcls),'\n','\n'));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constrained least squares l2-l1-TV (nonisotropic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(strcat(sprintf('Waiting for SUnSAL-TV ...'),'\n'));

% for sparsity regularization
lambda_ = 5*10^-4; % [0.5*10^-4 10^-3 5*10^-3 0.01 0.05 0.1 0.3 0.5 1];
% for spatial regularization
lambda_TV_ = 0.01; %[0.5*10^-4 10^-3 5*10^-3 0.01 0.05 0.1 0.3 0.5 1];

for i=1:length(lambda_)
    for j=1:length(lambda_TV_)
    
    lambda = lambda_(i);
    lambda_TV = lambda_TV_(j);
    
    t = clock;
    [X_hat_tv_ni,res,rmse_ni] = sunsal_tv(Dict,Y,'MU',0.05,'POSITIVITY','yes','ADDONE','no', ...
                               'LAMBDA_1',lambda,'LAMBDA_TV', lambda_TV, 'TV_TYPE','niso',...
                               'IM_SIZE',[75,75],'AL_ITERS',200, 'TRUE_X', XT,  'VERBOSE','no');
    SRE_tv_ni = 20*log10(norm(XT,'fro')/norm(X_hat_tv_ni-XT,'fro'));
    RMSE_tv_ni = Compute_RMSE(XT,X_hat_tv_ni);
    t_tv_ni = etime(clock,t);
    fprintf(strcat(sprintf('TV : (lambda - lambda_TV) = (%7.7f - %7.7f) ; (SRE - RMSE) = (%2.4f - %2.4f)',lambda, lambda_TV, SRE_tv_ni,RMSE_tv_ni),'\n','\n'));
    
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLUP_Lap GLUP with graph Laplacian regularization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(strcat(sprintf('Waiting for GLUP-Lap ...'),'\n'));

% Step 1 - Generate the affinity matrix
% ---------------------------------------
dmin = 0.3; % maximum SQUARE distance between pixels
flag = 0; % 1 if u want spatial information
wind = 3; % if flag =1, size of window to create averaged image 
t = clock;
[A, xy] = Affinity_Thresh(h,N,Y,dmin);
t_A = etime(clock,t);
Permuting_matrix = permute_grid_image(5,5,5); % For reordering pixels in image
Ap = Permuting_matrix'*A*Permuting_matrix; 
% figure; imshow(mat2gray(Ap)); axis on; 

% Step 2 - partition pixels into k clusters 
% ------------------------------------------
k = 10; % number of desired clusters
t = clock;
IDX = AMY(A,k); % use AMY method for Graph cut
t_clstr = etime(clock,t);
IDX_image = reshape(IDX,h,w);
% figure; imshow(mat2gray(IDX_image)); colormap cool; colorbar

% Step 3 - Solve GLUP_Lap for each cluster 
% -----------------------------------------
rho = 0.05; % penalty parameter for ADMM 
tol = 1e-5; % tol for stopping criteria
Nitermax = 200; % for ADMM

% for sparsity regularization
mu_ = 5*10^-4; % [0.5*10^-5 0.5*10^-4 10^-3 5*10^-3 0.01 0.05 0.1 0.3 0.5 1];
% for spatial regularization
mu_Lap_ = 0.5; % [0.5*10^-5 0.5*10^-4 10^-3 5*10^-3 0.01 0.05 0.1 0.3 0.5 1];

for i=1:length(mu_)
    for j=1:length(mu_Lap_)
        
        X_hat_GLUP_Lap = zeros(n,N); % memory allocation
        Z_hat_GLUP_Lap = zeros(n,N); % memory allocation
        Y_hat_GLUP_Lap = zeros(n,N); % memory allocation

        mu =  mu_(i); 
        mu_Lap =  mu_Lap_(j); 
        t = clock;
        
        for ii=1:k % loop for each cluster 
            %ii
            tmp = ((IDX-ii) == 0); % find pixels belonging to cluster ii
            indx_i = find(tmp);  clear tmp; % find pixels belonging to cluster ii
            Y_sub = Y(:,indx_i); % select spectra of those pixels
            A_sub = A(indx_i,indx_i); % select corresponding adjacency matrix
            XT_sub = XT(:,indx_i); % 
            X_hat_fcls_sub = X_hat_fcls(:,indx_i);
            [x_sub, z_sub, y_sub, RMSE_sub, t_sub, Niter ] = GLUP4_Lap(Dict,Y_sub,A_sub,rho,mu,mu_Lap,tol,Nitermax,XT_sub,X_hat_fcls_sub);
            X_hat_GLUP_Lap(:,indx_i) = x_sub;
            Z_hat_GLUP_Lap(:,indx_i) = z_sub;
            Y_hat_GLUP_Lap(:,indx_i) = y_sub;
        end
        
        t_glup_lAP = etime(clock,t);
        SRE_lapX = 20*log10(norm(X_hat_GLUP_Lap,'fro')/norm(X_hat_GLUP_Lap-XT,'fro'));
        SRE_lapZ = 20*log10(norm(Z_hat_GLUP_Lap,'fro')/norm(Z_hat_GLUP_Lap-XT,'fro'));
        SRE_lapY = 20*log10(norm(Y_hat_GLUP_Lap,'fro')/norm(Y_hat_GLUP_Lap-XT,'fro'));
        
        RMSE_lapX = Compute_RMSE(XT,X_hat_GLUP_Lap);
        
        fprintf(strcat(sprintf('Lap : (mu - mu_Lap) = (%7.7f - %7.7f) ; (SREX - RMSEX ) = (%7.3f = %2.4f)',mu, mu_Lap,SRE_lapX,RMSE_lapX),'\n','\n'));
    end
end

%%
%%%%%%%%%%%%%%%%%%%
%Displaying Figures
%%%%%%%%%%%%%%%%%%%

% Display the estimated abundance maps for each endmember
% Top to bottom:  true ones -> FCLS -> l2-l1 TV -> GLUP_Lap 
x_fcls_ = reshape(X_hat_fcls(1:P,:,:)',h,w,P); % reshape matrix into image
x_TV_ = reshape(X_hat_tv_ni(1:P,:,:)',h,w,P);
x_lap_ = reshape(X_hat_GLUP_Lap(1:P,:,:)',h,w,P);

figure; 
PP = 5;
for ic=1:PP
    subplot(4,PP,ic);imshow((x_fcls_(:,:,ic)));title('fcls')
    subplot(4,PP,ic+PP);imshow((x_TV_(:,:,ic)));title('tv')
    subplot(4,PP,ic+2*PP);imshow((x_lap_(:,:,ic)));title('lap')
    subplot(4,PP,ic+3*PP);imshow((image_X(:,:,ic)));title('true')
end
colormap pink;


PPP = P-PP;
if PPP~=0; figure; end
for ic=1:PPP
    subplot(4,PPP,ic);imshow((x_fcls_(:,:,ic+PP)));title('fcls')
    subplot(4,PPP,ic+PPP);imshow((x_TV_(:,:,ic+PP)));title('tv')
    subplot(4,PPP,ic+2*PPP);imshow((x_lap_(:,:,ic+PP)));title('lap')
    subplot(4,PPP,ic+3*PPP);imshow((image_X(:,:,ic+PP)));title('true')
end
colormap pink;

% % True image at Band 1
% Y_true = M*X; 
% Y_true_image = reshape(Y_true',h,w,L);
% figure; imshow((Y_true_image(:,:,1))); colormap pink ; colorbar
% 
% % Endmembers 
% figure; plot(Dict(:,supp)); axis([1 224 0 1])
% xlabel('Frequency band number','Interpreter','Latex','FontSize',16)
% ylabel('Reflectance','Interpreter','Latex','FontSize',16)

t_demo = etime(clock,t_demo); 
fprintf(strcat(sprintf('Total demo running time : %8.4f seconds',t_demo),'\n'));
