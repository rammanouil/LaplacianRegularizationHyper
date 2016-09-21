
function [X, Z, Y, RMSE, t, Niter ] = GLUP4_Lap(Sw,S,Adj,rho,mu,mu_L,tol,Nitermax,true_x,X0)

% x_sub, z_sub, y_sub, RMSE_sub, t_sub, Niter 
% ADMM algorithl for group LASSO with positivity constraint and sum-to-one constraint
%     variable X
%     minimize( square_pos(norm(S-SwX,'fro')) + lambda * sum(norms(X,2,2)) );
%     subject to
%         X >= nul;
%         uno * X == uno;

L = diag(sum(Adj)) - Adj;

[lSw,cSw] = size(Sw);
[lS,cS] = size(S);

if lS ~= lSw
    error('taille de Sw et S')
end

% init and precompute matrices

A = [eye(cSw,cSw);ones(1,cSw)];
B = [-eye(cSw,cSw);zeros(1,cSw)];
C = [zeros(cSw,cS);ones(1,cS)];

invSw2A2I = (Sw'*Sw + rho*(A'*A) + rho*eye(cSw))^-1;
SwS = Sw'*S;
AB = A'*B;
AC = A'*C;
invLrhoI = (2*mu_L*L + rho*eye(cS))^-1;

X = X0;
Y = X;
Z = X;

R = ones(cS+1,cS);

Lambda1 = zeros(cSw+1,cS);
Lambda2 = zeros(cSw,cS);

Niter = 0;

tol1 = sqrt(cS)*tol; 

t = clock; 
while norm(R) > tol1
    
    % X - minimization
    X = invSw2A2I*( SwS - A'*Lambda1 - rho*AB*Z + rho*AC - Lambda2 + rho*Y);
    
    % Z - minimization
    for k = 1:cSw;
        
        v = X(k,:) + (1/rho) * Lambda1(k,:);
        alpha = mu/rho;
        v_plus = max(v,0);
        normv_plus = norm(v_plus);
        
        if normv_plus < alpha
            Z(k,:) = zeros(1,cS);
        else
            Z(k,:) = (1 - alpha/normv_plus) * v_plus;
        end
        
    end
    
    % Y - minimization
    Y = (Lambda2 + rho*X)*invLrhoI;
    
    % Update Lagrange multipliers
    Lambda1 = Lambda1 + rho*(A*X + B*Z -C);
    Lambda2 = Lambda2 + rho*(X - Y);
    
    % Compute residual
    R = norm(S-Sw*X,'fro') + norm(X-Y,'fro') + norm(X-Z,'fro'); 
    
    Niter = Niter + 1;
    
    if Niter > Nitermax
        break 
    end 
   
    %disp(norm((X-true_x),'fro'));
    %disp(norm((Y-true_x),'fro'));
    %disp(norm((Z-true_x),'fro'));
end

RMSE_X = Compute_RMSE(X,true_x);
RMSE_Y = Compute_RMSE(Y,true_x);
RMSE_Z = Compute_RMSE(Z,true_x);

RMSE = min(RMSE_X,min(RMSE_Y,RMSE_Z));
t = etime(clock,t);
