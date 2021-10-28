% parameter selection for JSNMF
function [alpha, gamma, Inits] = parameter_selection(X1, X2,label)

num_clu = length(unique(label));
[W1,H1] = nndsvd(X1,num_clu,0); [W2,H2] = nndsvd(X2,num_clu,0);
Inits.W1 = W1; Inits.W2 = W2; Inits.H1 = H1; Inits.H2 = H2;

D1 = dist2(X1',X1'); S1 = affinityMatrix(D1,20); 
D2 = dist2(X2',X2'); S2 = affinityMatrix(D2,20); 

Inits.A1 = S1; Inits.D1 = diag(sum(S1)); Inits.L1 = Inits.D1-Inits.A1;
Inits.A2 = S2; Inits.D2 = diag(sum(S2)); Inits.L2 = Inits.D2-Inits.A2;

W = SNF({S1,S2},20);  % initialize S using SNF
Inits.S = W;
H1tH1 = H1'*H1; H2tH2 = H2'*H2;

err1 = norm(X1-W1*H1,'fro')^2; 
err2 = norm(X2-W2*H2,'fro')^2; 
err3 = 0.5*(norm(Inits.S-H1tH1,'fro')^2+norm(Inits.S-H2tH2,'fro')^2);
err4 = 0.5*trace(H1*Inits.L1*H1'+H2*Inits.L2*H2');
alpha = (err1+err2)/err3; 
gamma = abs((err1+err2)/err4); 

% scaling parameter
alpha = alpha/10;
gamma = gamma/10;

end