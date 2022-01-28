% parameter selection for JSNMF_batch_corr
function [alpha, gamma, Inits] = parameter_selection_batcheffect_pair5(X11, X12, X21, X23,X31,X34,X41, X45,X51,X56, num_factor1,num_factor2,num_factor3,num_factor4,num_factor5)
%num_clu = length(unique(label));
[W11,H11] = nndsvd(X11,num_factor1,0); [W2,H12] = nndsvd(X12,num_factor1,0);
[W21,H21] = nndsvd(X21,num_factor2,0); [W3,H23] = nndsvd(X23,num_factor2,0);
[W31,H31] = nndsvd(X31,num_factor3,0); [W4,H34] = nndsvd(X34,num_factor3,0);
[W41,H41] = nndsvd(X41,num_factor4,0); [W5,H45] = nndsvd(X45,num_factor4,0);
[W51,H51] = nndsvd(X51,num_factor5,0); [W6,H56] = nndsvd(X56,num_factor5,0);

Inits.W1 = (W11+W21+W31+W41+W51)/5; Inits.W2 = W2; Inits.W3= W3; Inits.W4 = W4; Inits.W5= W5; Inits.W6= W6;
Inits.H11 = H11; Inits.H12 = H12; Inits.H21 = H21; Inits.H23 = H23;
Inits.H31 = H31; Inits.H34 = H34; Inits.H41 = H41; Inits.H45 = H45;
Inits.H51 = H51; Inits.H56 = H56;

D11 = dist2(X11',X11'); S11 = affinityMatrix(D11,20); D12 = dist2(X12',X12'); S12 = affinityMatrix(D12,20); 
D21 = dist2(X21',X21'); S21 = affinityMatrix(D21,20); D23 = dist2(X23',X23'); S23 = affinityMatrix(D23,20); 
D31 = dist2(X31',X31'); S31 = affinityMatrix(D31,20); D34 = dist2(X34',X34'); S34 = affinityMatrix(D34,20); 
D41 = dist2(X41',X41'); S41 = affinityMatrix(D41,20); D45 = dist2(X45',X45'); S45 = affinityMatrix(D45,20); 
D51 = dist2(X51',X51'); S51 = affinityMatrix(D51,20); D56 = dist2(X56',X56'); S56 = affinityMatrix(D56,20); 

Inits.A11 = S11; Inits.D11 = diag(sum(S11)); Inits.L11 = Inits.D11-Inits.A11; Inits.A12 = S12; Inits.D12 = diag(sum(S12)); Inits.L12 = Inits.D12-Inits.A12;
Inits.A21 = S21; Inits.D21 = diag(sum(S21)); Inits.L21 = Inits.D21-Inits.A21; Inits.A23 = S23; Inits.D23 = diag(sum(S23)); Inits.L23 = Inits.D23-Inits.A23;
Inits.A31 = S31; Inits.D31 = diag(sum(S31)); Inits.L31 = Inits.D31-Inits.A31; Inits.A34 = S34; Inits.D34 = diag(sum(S34)); Inits.L34 = Inits.D34-Inits.A34;
Inits.A41 = S41; Inits.D41 = diag(sum(S41)); Inits.L41 = Inits.D41-Inits.A41; Inits.A45 = S45; Inits.D45 = diag(sum(S45)); Inits.L45 = Inits.D45-Inits.A45;
Inits.A51 = S51; Inits.D51 = diag(sum(S51)); Inits.L51 = Inits.D51-Inits.A51; Inits.A56 = S56; Inits.D56 = diag(sum(S56)); Inits.L56 = Inits.D56-Inits.A56;

snf_W1 = SNF({S11,S12},20); snf_W2 = SNF({S21,S23},20); 
snf_W3 = SNF({S31,S34},20); snf_W4 = SNF({S41,S45},20); snf_W5 = SNF({S51,S56},20);

Inits.S1 = snf_W1; Inits.S2 = snf_W2; Inits.S3 = snf_W3; Inits.S4 = snf_W4; Inits.S5 = snf_W5;
H11tH11 = H11'*H11; H12tH12 = H12'*H12;
H21tH21 = H21'*H21; H23tH23 = H23'*H23;
H31tH31 = H31'*H31; H34tH34 = H34'*H34;
H41tH41 = H41'*H41; H45tH45 = H45'*H45;
H51tH51 = H51'*H51; H56tH56 = H56'*H56;

err11 = norm(X11-W11*H11,'fro')^2; err12 = norm(X12-W2*H12,'fro')^2; 
err21 = norm(X21-W21*H21,'fro')^2; err23 = norm(X23-W3*H23,'fro')^2;
err31 = norm(X31-W31*H31,'fro')^2; err34 = norm(X34-W4*H34,'fro')^2; 
err41 = norm(X41-W41*H41,'fro')^2; err45 = norm(X45-W5*H45,'fro')^2; 
err51 = norm(X51-W51*H51,'fro')^2; err56 = norm(X56-W6*H56,'fro')^2; 

errS1 = 0.5*(norm(Inits.S1-H11tH11,'fro')^2+norm(Inits.S1-H12tH12,'fro')^2);
errS2 = 0.5*(norm(Inits.S2-H21tH21,'fro')^2+norm(Inits.S2-H23tH23,'fro')^2);
errS3 = 0.5*(norm(Inits.S3-H31tH31,'fro')^2+norm(Inits.S3-H34tH34,'fro')^2);
errS4 = 0.5*(norm(Inits.S4-H41tH41,'fro')^2+norm(Inits.S4-H45tH45,'fro')^2);
errS5 = 0.5*(norm(Inits.S5-H51tH51,'fro')^2+norm(Inits.S5-H56tH56,'fro')^2);

err7 = 0.25*trace(H11*Inits.L11*H11'+H12*Inits.L12*H12'+H21*Inits.L21*H21'+H23*Inits.L23*H23');
err8 = 0.25*trace(H31*Inits.L31*H31'+H34*Inits.L34*H34'+H41*Inits.L41*H41'+H45*Inits.L45*H45');
err9 = 0.25*trace(H51*Inits.L51*H51'+H56*Inits.L56*H56');

tot1 = err11+err12+err21+err23+err31+err34+err41+err45+err51+err56;
alpha = tot1/(errS1+errS2+errS3+errS4+errS5); 
gamma = abs(tot1/(err7+err8+err9)); 

% scaling parameter
alpha = alpha/10;
gamma = gamma/10;

end