function [W1,W2,W3,H1,H2,H3,S,iter,objs] = jsnmf_3mod(X1,X2,X3,alpha,gamma,Inits)
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


Maxiter = 500; 
W1 = Inits.W1; W2 = Inits.W2; W3 = Inits.W3; H3 = Inits.H3'; H1 = Inits.H1'; H2 = Inits.H2';
A1 = Inits.A1; D1 = Inits.D1; S = Inits.S; S = S./repmat(sum(S),size(S,1),1);
A2 = Inits.A2; D2 = Inits.D2; A3 = Inits.A3; D3 = Inits.D3; n = size(S,1); 

theta1 = 1/3; theta2 = 1/3; theta3 = 1/3;
obj_old = 1; stop_rule = 2; yita1 = 1; yita2 = 0.5; yita3 = 0.5; lamda = 0.5; %yita1: 1; yita2: 0.5;
objs = zeros(Maxiter,1);
clear Inits

disp('iteration starts!');
for iter = 1:Maxiter
    % update W1 W2 using bpp algorithm defaulted, faster when the data volumne is large
    [W1,~,~,~] = nnlsm_blockpivot(H1, X1', 0, W1');
    [W2,~,~,~] = nnlsm_blockpivot(H2, X2', 0, W2');
    [W3,~,~,~] = nnlsm_blockpivot(H3, X3', 0, W3');
    W1 = W1'; W2 = W2'; W3 = W3';
    
    % updating w1 w2 via multiplication rule
    % H1tH1 = H1'*H1; H2tH2 = H2'*H2;
    % W1 = W1.*((X1*H1)./max(W1*H1tH1,1e-10));
    % W2 = W2.*((X2*H2)./max(W2*H2tH2,1e-10));
    
    % update H1 H2
    W1tW1 = W1'*W1; W2tW2 = W2'*W2; W3tW3 = W3'*W3; H1tH1 = H1'*H1; H2tH2 = H2'*H2; H3tH3 = H3'*H3; 
    H1 = H1.*((alpha*S'*H1+X1'*W1+gamma*theta1.^2*A1*H1+yita1*H1)./max(gamma*theta1.^2*D1*H1+(alpha+yita1)*H1*H1tH1+H1*W1tW1,1e-10));
    H2 = H2.*((alpha*S'*H2+X2'*W2+gamma*theta2.^2*A2*H2+yita2*H2)./max(gamma*theta2.^2*D2*H2+(alpha+yita2)*H2*H2tH2+H2*W2tW2,1e-10));
    H3 = H3.*((alpha*S'*H3+X3'*W3+gamma*theta3.^2*A3*H3+yita3*H3)./max(gamma*theta3.^2*D3*H3+(alpha+yita3)*H3*H3tH3+H3*W3tW3,1e-10));

    % update S
    H1H1t = H1*H1'; H2H2t = H2*H2'; H3H3t = H3*H3'; Q = (H1H1t+H2H2t+H3H3t)/3;
    tot = repmat(sum(S),size(S,1),1);
    S = S.*((alpha*Q+lamda*ones(n))./max(alpha*S+lamda*tot,1e-10));    
    clear H1H1t H2H2t Q H1tH1 H2tH2 X2SR W1tW1 W2tW2 H3H3t
   
    % updata theta
    tmp1 = 1/trace(H1'*(D1-A1)*H1); tmp2 = 1/trace(H2'*(D2-A2)*H2); tmp3 = 1/trace(H3'*(D3-A3)*H3);
    tmp = tmp1+tmp2+tmp3;
    theta1 = tmp1/tmp; theta2 = tmp2/tmp; theta3 = tmp3/tmp;
    
    if stop_rule == 2
        obj = compute_obj(X1,X2,X3,S,W1,W2,W3,H1,H2,H3,D1,A1,D2,A2,D3,A3,alpha,gamma,theta1,theta2,theta3);
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
S = (S+S')/2;
end

function obj = compute_obj(X1,X2,X3,S,W1,W2,W3,H1,H2,H3,D1,A1,D2,A2,D3,A3,alpha,gamma,theta1,theta2,theta3)
  L1 = D1-A1; L2 = D2-A2; L3 = D3-A3;
  obj = norm(X1-W1*H1','fro')^2+norm(X2-W2*H2','fro')^2+norm(X3-W3*H3','fro')^2+0.5*alpha*(norm(S-H1*H1','fro')^2+...
     norm(S-H2*H2','fro')^2+norm(S-H3*H3','fro')^2)+gamma*(trace(theta1.^2*H1'*L1*H1+theta2.^2*H2'*L2*H2+theta3.^2*H3'*L3*H3));
end

