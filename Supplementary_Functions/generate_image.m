function [r,M,a,std_noise]=generate_image(M_database, R, model, N, SNR, a)
% Image Generation for Test
% ----- Input parameters ----- 
% M_database: database selection
% R: endmember 
% model: linear, bilinear, intimate, postnonlinear
% N: number of pixels
% SNR: signal-to-noise ratio
% a: input abundances.(if provided)
% ---- Outputs ------
% r: the mixed image
% a: abundance used


if nargin < 6;
   % use genearated abundances
end

% endmember database choice
if M_database == 1
    load ../material.mat
    M = [allunite(:,2) budd(:,2) calcite(:,2) kaolinite(:,2) muscovite(:,2)];
    M = M(:,1:R);
elseif M_database == 2
    load M_CD_536_6.mat
    M = M(:,1:R);
    M = M/100;
elseif M_database == 3 
    load endmembers.mat
    M = M(:,2:end);
    M = M(:,1:R);
    % M= M(:,[3,5,6]);
    % M(:,end)=M(:,end)*0.1;
elseif M_database == 4
    load Albedo.mat
    M = Albedo(:,[1,3,4]);
elseif M_database == 5
    load ENVI_Min_420_8.mat
    % M = M(:,2:end);
    % M = M(:,1:R);
    M = M(:,[3,5,6]);
elseif M_database == 100
    %M = eye(R);
    M=[0.9, 0.6, 0.3;
       0.8, 0.9, 0.2;
       0.5, 0.7, 0.8];
elseif M_database == 33 
    load endmembersTV_article.mat
    M = M(:,1:R);
    % M= M(:,[3,5,6]);
    % M(:,end)=M(:,end)*0.1;
elseif M_database == 44 
    load M9.mat
    M = M(:,1:R);
    % M= M(:,[3,5,6]);
    % M(:,end)=M(:,end)*0.1;
else
    %for DC2
    load DATA_DC2;
    M = A(:,selected);
  %  M=A;
     M = M(:,1:R);
end
%M(:,2)=-M(:,2);
%M = M(1:2:end,:);
%M = chooseM(M,0.9995);



[L,R] = size(M);
%figure, plot(M);

if nargin < 6
    % Abundance loading
      a=AbundanceGen(R,N);
end

   % a=AbundanceGen(R,2500,0.8);
% Image generation
[r_pure,a] = hypermix(M,N,model,a);

% Additive noise
pw_signal = norm(r_pure,'fro')^2/(N*L);
pw_noise = pw_signal/(10^(SNR/10));
std_noise = sqrt(pw_noise);
noise = std_noise*randn(L,N);
%
r = r_pure + noise;
end

