% Code for reprodcucing Figure 3 
% from the manuscript Malioutov & Slavov, Convex Total Least Squares 
% http://proceedings.mlr.press/v32/malioutov14.html


addpath slra-0.5
%%
n=120;
sym_dat = -5*ones( n, 3*4 ); 

Errors = [3 5]; %7.8
clear Corrs s

for Noise = 1:2
    Error = Errors(Noise);
    col = (Noise-1)*8;
    
    for Iter = 1:n
            P = [1 2 3
                 5 3 1];

            Z = diag( [1 5 8] );
            S = [1 0
                 0 1
                 1 1];


             % The noisless data should be
             X = Z * S * P;

             % Add noise to the data
             Xn = X;
             Xn(1) = Error; 
             %Xn(end) = Xn(end) - 2*Error;
             %Xn(3,1) = Xn(end) + 3*Error;
             %Xn(3,2) = Xn(3,2) - 3*Error;
             Outliers_Support = X ~= Xn;
             X = Xn + 0.1*randn( size(X) );



            [M, N] = size(X);
            [~, K] = size(S);

            A1 = kron(S, eye(N));
            A2 = zeros(M*N, M);
            for i = 1:M
               A2((i-1)*N + (1:N), i) = -X(i, :)';
            end
            A = [A1 A2]; 


            % 1 Solve by SVD
            [u, ss, v] = svd( A ); 
            col_dim = sum( diag(ss) >1e-10 );
            num_dim = size(A,2)-col_dim; 

            hat_P = zeros(size(P));
            for i = 1:K
                hat_P(i, :) = v( (i-1)*N + (1:N), end )';
            end
            if  median( sign( hat_P(:) ) ) == 0
                hat_P = hat_P/hat_P(1);
            else
                hat_P = hat_P * median( sign( hat_P(:) ) ); 
            end

            Corrs.svd = diag( corr( P', hat_P' ) );
            
            
             %2 EM SVDimpute
               %EM Loop
               X_imputted = X;
               for jj=1:1e3
                   
                    [u, ss, v] = svd( X_imputted );
                    
                    %vt = v(1:end-1,:)';
                    ut = u(:,1:end-1);
                    sol = linsolve( ut, X_imputted );
                   
                    X_imputted_new = ut*sol;

                    
                    if  norm( X_imputted_new(Outliers_Support) - X_imputted(Outliers_Support), 'fro' ) < 1e-5, break, end
                    X_imputted(Outliers_Support) = X_imputted_new(Outliers_Support);
               end
                for i = 1:M
                    A2_imputted((i-1)*N + (1:N), i) = -X_imputted(i, :)';
                    %Weights((i-1)*N + (1:N), i) = Outliers_Support(i, :)'==0;
                end
                A_imputted = [A1 A2_imputted];
                %W = [ones( size(A1) ) Weights];
                
             
             
            [~, ss, v] = svd( A_imputted ); 
            col_dim = sum( diag(ss) >1e-10 );
            num_dim = size(A,2)-col_dim; 

            hat_P = zeros(size(P));
            for i = 1:K
                hat_P(i, :) = v( (i-1)*N + (1:N), end )';
            end
            if  median( sign( hat_P(:) ) ) == 0
                hat_P = hat_P/hat_P(1);
            else
                hat_P = hat_P * median( sign( hat_P(:) ) ); 
            end
 

            Corrs.SVDimpute = diag( corr( P', hat_P' ) );
            


            % 3 Solve by slra-0.5 and SVD
            p = A; %rand( 5 );
            s.m = ones( size(p,1), 1);
            %S.w = W; 
            s.w = ones( size(p) );
            s.w(1,end-2) = 0;
            %s.w(end,end) = 0;
               %s.w(end-2,end) = 0;
            %s.w(end-1,end) = 0;   
            r = 8;

            [ph, info] = slra(p, s, r );
            ph = reshape( ph, size(p) );

            [u, ss, v] = svd( ph ); 
            col_dim = sum( diag(ss) >1e-10 );
            num_dim = size(A,2)-col_dim; 

            hat_P = zeros(size(P));
            for i = 1:K
                hat_P(i, :) = v( (i-1)*N + (1:N), end )';
            end
            if  median( sign( hat_P(:) ) ) == 0
                hat_P = hat_P/hat_P(1);
            else
                hat_P = hat_P * median( sign( hat_P(:) ) ); 
            end

            Corrs.slra = diag( corr( P', hat_P' ) );


            % 4 Solve by convex STLS RWNN
            if 0 || n==1
                [x_rec3, X_tilde] = RWNN_tls(A, s.w );

                [u, ss, v] = svd( X_tilde ); 
                col_dim = sum( diag(ss) >1e-10 );
                num_dim = size(A,2)-col_dim; 

                hat_P = zeros(size(P));
                for i = 1:K
                    hat_P(i, :) = v( (i-1)*N + (1:N), end )';
                end
                if  median( sign( hat_P(:) ) ) == 0
                    hat_P = hat_P/hat_P(1);
                else
                    hat_P = hat_P * median( sign( hat_P(:) ) ); 
                end

                Corrs.RWNN = diag( corr( P', hat_P' ) );
            else
                Corrs.RWNN = 1;
            end



           sym_dat(Iter,1+col) =  mean( Corrs.svd );
           sym_dat(Iter,2+col) =  mean( Corrs.SVDimpute );
           sym_dat(Iter,3+col) =  mean( Corrs.slra );
           sym_dat(Iter,4+col) =  mean( Corrs.RWNN );

    end
end





%
if size(sym_dat,1) ==1, sym_dat=[sym_dat;sym_dat]; end

boxplot( sym_dat,  'color', 'gbkr', 'notch', 'off', 'whisker', 1.5 )%

ylim( [-1.1 1.1] )
set( gcf, 'Units', 'normalized' )
set( gcf, 'Position', [0.045    0.348    0.521    0.51] )
set( gca, 'Position', [0.13  0.145 0.84 0.83] )
set( gca, 'FontSize', 14, 'FontWeight', 'Bold' );

set( gca, 'Xtick', [2 8],  'Xticklabel', { 'Low Noise', 'High Noise' } );



h(1) = ylabel( 'Inference Accuracy, $\langle U, \hat U \rangle$' );
h(2) = text( 4.5,  0.3, 'SVD', 'color', 'g' );
h(3) = text( 4.5,  0.1,   'SVDimpute EM', 'color', 'b'  );
h(4) = text( 4.5, -0.1,   'Non-Convex STLS', 'color', 'k'  );
h(5) = text( 4.5, -0.3, 'Convex STLS', 'color', 'r' );
set( h,  'fontsize',       22,...
         'fontweight',    'bold' );
set( h,  'interpreter',    'latex'  ); 

%%
pdf( 'Suggested_Figure', [], 1 )



















