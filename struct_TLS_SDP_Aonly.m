function A_tilde = struct_TLS_SDP_Aonly(A, Msk, max_iter, W_e)

%%% find a closest rank-deficient matrix A_tilde to A, 
%%% where we're only allowed to perturb entries in Msk.
%%% similar to STLS -- but without the rhs -- we're looking
%%% for nullspace elements of structured perturbations of A.
%%% W_e is entry-wise weight on the frobenius norm

addpaths;

if nargin < 3,  max_iter = 2; end %%% reweighting iterations
if nargin < 4, W_e = ones(size(A)); end
%%% if max_iter = 1 -- we have non-reweighted TLS...

for iter = 1:max_iter
    fprintf('------------------->\n reweighting iter %d\n -------------------->\n', iter);
    if iter == 1
        [dA, lst_beta, lst_Y, lst_Z] = solve_STLS(A, Msk, W_e);
        %[E, lst_beta_tilde, Wy_new, Wz_new] = solve_STLS_ALM(A, Msk);
    else
        [dA, lst_beta, lst_Y, lst_Z] = solve_STLS(A, Msk, W_e, lst_beta, lst_Y, lst_Z);
        %[E, lst_beta_tilde, Wy_new, Wz_new] = solve_STLS_ALM(A, Msk, lst_beta_tilde, Wy_new, Wz_new);
    end

    A_tilde = A - dA;
    
end
    
disp('');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dX_rec, beta, Y_last, Z_last] = solve_STLS(X, Msk, W_e, beta, Y, Z)
%%% solve an STLS problem by reweighting
%%%
%%% this function allows reweighted-trace-norm (i.e. resolving from the previous
%%% solution -- if we're given Y, Z, beta_old). It returns the new Y, Z,
%%% beta

[M, N] = size(X);
if nargin <= 3, %%% reweighting starts only after the first iteration
    % Y = eye(M); Z = eye(N);   
    WW = eye(M+N);  
    beta = 1 ; %beta = 0.9;
else
    delta = 0.01;
    WW = inv(blkdiag(Y, Z) + delta*eye(M+N));
end

inds_err = find(Msk);

options = sdpsettings('solver', 'SDPT3', 'verbose', 0);

X_sdp = sdpvar(M, N, 'full');
Y_sdp = sdpvar(M);
Z_sdp = sdpvar(N);

XX = [Y_sdp (X-X_sdp); (X-X_sdp)' Z_sdp];
F = [XX >= 0, X_sdp(inds_err) == 0];

iter = 1; max_iter = 100; beta_mn = 0; beta_mx = Inf;
last_good_beta = Inf; dbeta = Inf;

zr_thresh = 5e-6;  %%% judge if sv is zero or not!

while dbeta > min(0.01, 0.25*beta) && iter < max_iter
    %    obj2 = 1/2*beta*(norm(WW*XX, 'nuclear')) + norm(W_e.*X_sdp, 'fro')^2;
    obj2 = 1/2*beta*(norm(WW*XX, 'nuclear')) + sum(sum((W_e.*X_sdp).^2));
    status = solvesdp(F, obj2, options);
    dX_rec = double(X_sdp);    
    X_lr_sdp = X - dX_rec;

    fprintf('iter %d: beta = %.5f, beta_rg: [ %.5f, %.5f], svd of X-X_rec: \n', ...
        iter, beta, beta_mn, beta_mx);
    s = svd(X_lr_sdp);
    %disp(s);    
    s_min = min(s);
        
    if s_min > zr_thresh
        beta_mn = beta; beta = 2*beta; 
        if ~isinf(beta_mx), beta = 0.5*(beta_mn + beta_mx); end
    end
        
    if s_min <= zr_thresh  %%% valid solution
        last_good_beta = beta; %%% remember last valid solution (best seen so far)
        lg_X = dX_rec; lg_Y = double(Y_sdp); lg_Z = double(Z_sdp);
        
        beta_mx = beta; beta = 0.5*beta;
        if ~isinf(beta_mn), beta = 0.5*(beta_mn + beta_mx); end
    end
    
    dbeta = beta_mx - beta_mn;  %%% beta-gap
    
    iter = iter + 1;
end

fprintf('last-good-beta: %.5f, beta-rg: [ %.5f, %.5f]\n', last_good_beta, beta_mn, beta_mx);

beta = last_good_beta;
dX_rec = lg_X; Y_last = lg_Y; Z_last = lg_Z;

disp('');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addpaths

addpath(genpath('Yalmip'));
addpath(genpath('SDPT3-4.0'));
addpath(genpath('slra-0.5'));
