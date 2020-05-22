function [x_rec3, X_tilde] = RWNN_tls(X, s)

%%% unweighted RWNN-TLS


Msk = zeros(size(X)); 
max_iter = 4;

X_tilde = struct_TLS_SDP_Aonly(X, Msk, max_iter, s);
[U, S, V] = svd(X_tilde);
n_rec3 = V(:, end);
n_rec3 = n_rec3/n_rec3(end);
x_rec3 = n_rec3(1:(end-1));

disp('');
