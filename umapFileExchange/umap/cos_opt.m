function W = cos_opt(X)
%%
W = X*X';
DX = sqrt(diag(W))*sqrt(diag(W))';
W = W./DX;

end
