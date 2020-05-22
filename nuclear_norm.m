function y = nuclear_norm ( x )

[~, s, ~] = svd( x );
y=sum(diag(s));