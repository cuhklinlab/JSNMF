function [W1,W2,W3,W4,W5,W6,H11,H12,H21,H23,H31,H34,H41,H45,H51,H56,S1,S2,S3,S4,S5,iter,objs] = jsnmf_batch_correct_pari5(X11,X12,X21,X23,X31,X34,X41, X45,X51,X56,alpha,gamma,Inits)
% the function is used to implement JSNMF algorithm
% input:
% X1: p*n normalized single-cell transcriptomic data, row denotes gene and cloumn denotes cell
% X2: q*n normalized single-cell epigenomic data, row denotes region and column denotes cell
% alpha, gamma: parameter
% Inits: initializations of W1, W2, H1, H2, A1, A2 and S
% output:
% W1, W2, H1, H2: the low-dimensional factor loading for genens, regions, cells(H1 for RNA, H2 for epigenomic)
% S: the cell-cell similarity matrix, used to cluster cell subpopulations, visulize and so on
% iter: the number of iteration
% objs: the value of objective functions when it converges


Maxiter = 2000; 
W1 = Inits.W1; W2 = Inits.W2; W3 = Inits.W3; W4 = Inits.W4; W5 = Inits.W5; W6 = Inits.W6;
H11 = Inits.H11'; H12 = Inits.H12'; H21 = Inits.H21'; H23 = Inits.H23';
A11 = Inits.A11; D11 = Inits.D11;  A12 = Inits.A12; D12 = Inits.D12; 
A21 = Inits.A21; D21 = Inits.D21;  A23 = Inits.A23; D23 = Inits.D23;
S1 = Inits.S1; S1 = S1./repmat(sum(S1),size(S1,1),1); n1 = size(S1,1); 
S2 = Inits.S2; S2 = S2./repmat(sum(S2),size(S2,1),1); n2 = size(S2,1); 

H31 = Inits.H31'; H34 = Inits.H34'; H41 = Inits.H41'; H45 = Inits.H45';
A31 = Inits.A31; D31 = Inits.D31;  A34 = Inits.A34; D34 = Inits.D34; 
A41 = Inits.A41; D41 = Inits.D41;  A45 = Inits.A45; D45 = Inits.D45;
S3 = Inits.S3; S3 = S3./repmat(sum(S3),size(S3,1),1); n3 = size(S3,1); 
S4 = Inits.S4; S4 = S4./repmat(sum(S4),size(S4,1),1); n4 = size(S4,1); 

H51 = Inits.H51'; A51 = Inits.A51; D51 = Inits.D51;
H56 = Inits.H56'; A56 = Inits.A56; D56 = Inits.D56;
S5 = Inits.S5; S5 = S5./repmat(sum(S5),size(S5,1),1); n5 = size(S5,1); 

theta11 = 1/10; theta12 = 1/10; theta21 = 1/10; theta23 = 1/10;
theta31 = 1/10; theta34 = 1/10; theta41 = 1/10; theta45 = 1/10;
theta51 = 1/10; theta56 = 1/10;

obj_old = 1; stop_rule = 2; yita11 = 1; yita12 = 0.5; yita21 = 1; yita23 = 0.5; lamda = 0.5; %yita1: 1; yita2: 0.5;
yita31 = 1; yita34 = 0.5; yita41 = 1; yita45 = 0.5; yita51 = 1; yita56 = 0.5;
objs = zeros(Maxiter,1);
clear Inits

disp('iteration starts!');
for iter = 1:Maxiter
    % update W1 W2 using bpp algorithm defaulted, faster when the data volumne is large
    %[W1,~,~,~] = nnlsm_blockpivot(H1, X1', 0, W1');
    H11tH11 = H11'*H11; H21tH21 = H21'*H21; H31tH31 = H31'*H31; H41tH41 = H41'*H41; H51tH51 = H51'*H51;
    W1 = W1.*((X11*H11+X21*H21+X31*H31+X41*H41+X51*H51)./max(W1*(H11tH11+H21tH21+H31tH31+H41tH41+H51tH51),1e-10));
    [W2,~,~,~] = nnlsm_blockpivot(H12, X12', 0, W2');
    [W3,~,~,~] = nnlsm_blockpivot(H23, X23', 0, W3');
    [W4,~,~,~] = nnlsm_blockpivot(H34, X34', 0, W4');
    [W5,~,~,~] = nnlsm_blockpivot(H45, X45', 0, W5');
    [W6,~,~,~] = nnlsm_blockpivot(H56, X56', 0, W6');
    W2 = W2'; W3 = W3'; W4 = W4'; W5 = W5'; W6 = W6';
    
    % updating w1 w2 via multiplication rule
    % H1tH1 = H1'*H1; H2tH2 = H2'*H2;
    % W1 = W1.*((X1*H1)./max(W1*H1tH1,1e-10));
    % W2 = W2.*((X2*H2)./max(W2*H2tH2,1e-10));
    
    % update H1 H2
    W1tW1 = W1'*W1; W2tW2 = W2'*W2; W3tW3 = W3'*W3; H12tH12 = H12'*H12; H23tH23 = H23'*H23; 
    H11 = H11.*((alpha*S1'*H11+X11'*W1+gamma*theta11.^2*A11*H11+yita11*H11)./max(gamma*theta11.^2*D11*H11+(alpha+yita11)*H11*H11tH11+H11*W1tW1,1e-10));
    H12 = H12.*((alpha*S1'*H12+X12'*W2+gamma*theta12.^2*A12*H12+yita12*H12)./max(gamma*theta12.^2*D12*H12+(alpha+yita12)*H12*H12tH12+H12*W2tW2,1e-10));
    H21 = H21.*((alpha*S2'*H21+X21'*W1+gamma*theta21.^2*A21*H21+yita21*H21)./max(gamma*theta21.^2*D21*H21+(alpha+yita21)*H21*H21tH21+H21*W1tW1,1e-10));
    H23 = H23.*((alpha*S2'*H23+X23'*W3+gamma*theta23.^2*A23*H23+yita23*H23)./max(gamma*theta23.^2*D23*H23+(alpha+yita23)*H23*H23tH23+H23*W3tW3,1e-10));
    
    W4tW4 = W4'*W4; W5tW5 = W5'*W5; W6tW6 = W6'*W6; H34tH34 = H34'*H34; H45tH45 = H45'*H45; H56tH56 = H56'*H56;
    H31 = H31.*((alpha*S3'*H31+X31'*W1+gamma*theta31.^2*A31*H31+yita31*H31)./max(gamma*theta31.^2*D31*H31+(alpha+yita31)*H31*H31tH31+H31*W1tW1,1e-10));
    H34 = H34.*((alpha*S3'*H34+X34'*W4+gamma*theta34.^2*A34*H34+yita34*H34)./max(gamma*theta34.^2*D34*H34+(alpha+yita34)*H34*H34tH34+H34*W4tW4,1e-10));
    H41 = H41.*((alpha*S4'*H41+X41'*W1+gamma*theta41.^2*A41*H41+yita41*H41)./max(gamma*theta41.^2*D41*H41+(alpha+yita41)*H41*H41tH41+H41*W1tW1,1e-10));
    H45 = H45.*((alpha*S4'*H45+X45'*W5+gamma*theta45.^2*A45*H45+yita45*H45)./max(gamma*theta45.^2*D45*H45+(alpha+yita45)*H45*H45tH45+H45*W5tW5,1e-10));
    
    H51 = H51.*((alpha*S5'*H51+X51'*W1+gamma*theta51.^2*A51*H51+yita51*H51)./max(gamma*theta51.^2*D51*H51+(alpha+yita51)*H51*H51tH51+H51*W1tW1,1e-10));
    H56 = H56.*((alpha*S5'*H56+X56'*W6+gamma*theta56.^2*A56*H56+yita56*H56)./max(gamma*theta56.^2*D56*H56+(alpha+yita56)*H56*H56tH56+H56*W6tW6,1e-10));
        
    % update S
    H11H11t = H11*H11'; H12H12t = H12*H12'; H21H21t = H21*H21'; H23H23t = H23*H23';
    Q1 = 0.5*(H11H11t+H12H12t); Q2 = 0.5*(H21H21t+H23H23t);
    H31H31t = H31*H31'; H34H34t = H34*H34'; H41H41t = H41*H41'; H45H45t = H45*H45'; H51H51t = H51*H51'; H56H56t = H56*H56';
    Q3 = 0.5*(H31H31t+H34H34t); Q4 = 0.5*(H41H41t+H45H45t);
    Q5 = 0.5*(H51H51t+H56H56t);
    
    tot1 = repmat(sum(S1),size(S1,1),1); tot2 = repmat(sum(S2),size(S2,1),1);
    tot3 = repmat(sum(S3),size(S3,1),1); tot4 = repmat(sum(S4),size(S4,1),1); tot5 = repmat(sum(S5),size(S5,1),1); 
    
    S1 = S1.*((alpha*Q1+lamda*ones(n1))./max(alpha*S1+lamda*tot1,1e-10)); 
    S2 = S2.*((alpha*Q2+lamda*ones(n2))./max(alpha*S2+lamda*tot2,1e-10));
    S3 = S3.*((alpha*Q3+lamda*ones(n3))./max(alpha*S3+lamda*tot3,1e-10)); 
    S4 = S4.*((alpha*Q4+lamda*ones(n4))./max(alpha*S4+lamda*tot4,1e-10)); 
    S5 = S5.*((alpha*Q5+lamda*ones(n5))./max(alpha*S5+lamda*tot5,1e-10)); 
    clear H11H11t H12H12t Q1 Q2 H21H21t H23H23t H11tH11 H21tH21 H31H31t H34H34t H41H41t H45H45t H51H51t H56H56t Q3 Q4 Q5
   
    % updata theta
    tmp11 = 1/trace(H11'*(D11-A11)*H11); tmp12 = 1/trace(H12'*(D12-A12)*H12); 
    tmp21 = 1/trace(H21'*(D21-A21)*H21); tmp23 = 1/trace(H23'*(D23-A23)*H23);
    
    tmp31 = 1/trace(H31'*(D31-A31)*H31); tmp34 = 1/trace(H34'*(D34-A34)*H34); 
    tmp41 = 1/trace(H41'*(D41-A41)*H41); tmp45 = 1/trace(H45'*(D45-A45)*H45);
    tmp51 = 1/trace(H51'*(D51-A51)*H51); tmp56 = 1/trace(H56'*(D56-A56)*H56);

    tmp = (tmp11+tmp12)+(tmp21+tmp23)+(tmp31+tmp34)+(tmp41+tmp45)+(tmp51+tmp56);
    theta11 = tmp11/tmp; theta12 = tmp12/tmp;theta21 = tmp21/tmp; theta23 = tmp23/tmp;
    theta31 = tmp31/tmp; theta34 = tmp34/tmp;theta41 = tmp41/tmp; theta45 = tmp45/tmp;
    theta51 = tmp51/tmp; theta56 = tmp56/tmp;
    
    if stop_rule == 2
        obj = compute_obj(X11,X12,X21,X23,X31,X34,X41,X45,X51,X56,S1,S2,S3,S4,S5,W1,W2,W3,W4,W5,W6,H11,H12,H21,H23,H31,H34,H41,H45,H51,H56,D11,A11,D12,A12,D21,A21,D23,A23,D31,A31,D34,A34,D41,A41,D45,A45,D51,A51,D56,A56, alpha,gamma,theta11,theta12,theta21,theta23,theta31,theta34,theta41,theta45,theta51,theta56);
        objs(iter,1) = obj;
        if (abs(obj_old-obj)/obj_old < 10^(-5) && iter > 1) || iter == Maxiter
            disp('converged!');
            break;
        end
        obj_old = obj;
    end
    if mod(iter,10) == 0
        disp(['number of iteration:',num2str(iter),'  obj:',num2str(obj)]);
        %disp(obj);
    end
end
S1 = (S1+S1')/2; S2 = (S2+S2')/2; S3 = (S3+S3')/2; S4 = (S4+S4')/2; S5 = (S5+S5')/2;
end

function obj = compute_obj(X11,X12,X21,X23,X31,X34,X41,X45,X51,X56,S1,S2,S3,S4,S5,W1,W2,W3,W4,W5,W6,H11,H12,H21,H23,H31,H34,H41,H45,H51,H56,D11,A11,D12,A12,D21,A21,D23,A23,D31,A31,D34,A34,D41,A41,D45,A45,D51,A51,D56,A56, alpha,gamma,theta11,theta12,theta21,theta23,theta31,theta34,theta41,theta45,theta51,theta56)
  L11 = D11-A11; L12 = D12-A12;L21 = D21-A21; L23 = D23-A23;
  L31 = D31-A31; L34 = D34-A34;L41 = D41-A41; L45 = D45-A45;
  L51 = D51-A51; L56 = D56-A56;
  obj = norm(X11-W1*H11','fro')^2+norm(X12-W2*H12','fro')^2+norm(X21-W1*H21','fro')^2+norm(X23-W3*H23','fro')^2+...
        norm(X31-W1*H31','fro')^2+norm(X34-W4*H34','fro')^2+norm(X41-W1*H41','fro')^2+norm(X45-W5*H45','fro')^2+norm(X51-W1*H51','fro')^2+norm(X56-W6*H56','fro')^2+...
   0.5*alpha*(norm(S1-H11*H11','fro')^2+norm(S1-H12*H12','fro')^2+norm(S2-H21*H21','fro')^2+norm(S2-H23*H23','fro')^2+...
   norm(S3-H31*H31','fro')^2+norm(S3-H34*H34','fro')^2+norm(S4-H41*H41','fro')^2+norm(S4-H45*H45','fro')^2+norm(S5-H51*H51','fro')^2+norm(S5-H56*H56','fro')^2)+...
gamma*(trace(theta11.^2*H11'*L11*H11+theta12.^2*H12'*L12*H12+theta21.^2*H21'*L21*H21+theta23.^2*H23'*L23*H23)+...
trace(theta31.^2*H31'*L31*H31+theta34.^2*H34'*L34*H34+theta41.^2*H41'*L41*H41+theta45.^2*H45'*L45*H45+ theta51.^2*H51'*L51*H51+theta56.^2*H56'*L56*H56));
end

