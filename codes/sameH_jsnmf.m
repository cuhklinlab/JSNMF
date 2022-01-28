function [W1,W2,H,iter,objs] = sameH_jsnmf(X1,X2,gamma,Inits)
% H shared by scRNA and scATAC data
Maxiter = 300; 
W1 = Inits.W1; W2 = Inits.W2; H = Inits.H1'; 
A1 = Inits.A1; D1 = Inits.D1;
A2 = Inits.A2; D2 = Inits.D2; 
clear Inits

theta1 = 1/2; theta2 = 1/2;
obj_old = 1; stop_rule = 2; yita = 1; 
objs = zeros(Maxiter,1);
for iter = 1:Maxiter
    % update W1 W2
    [W1,~,~,~] = nnlsm_blockpivot(H, X1', 0, W1');
    [W2,~,~,~] = nnlsm_blockpivot(H, X2', 0, W2');
    W1 = W1'; W2 = W2';
   
    %updating w1 w2 via multiplication rule
%     HtH = H'*H;
%     W1 = W1.*((X1*H)./max(W1*HtH,1e-10));
%     SR = S; SR(Index) = 0; 
%     X2SR = X2*SR;
%     W2 = W2.*((X2*H)./max(W2*HtH,1e-10));
    
    % update H1 H2
    W1tW1 = W1'*W1; W2tW2 = W2'*W2; HtH = H'*H; 
    H = H.*(X1'*W1+X2'*W2+yita*H+gamma*(theta1.^2*A1+theta2.^2*A2)*H)./max(H*(W1tW1+W2tW2)+yita*H*HtH+gamma*(theta1.^2*D1+theta2.^2*D2)*H,1e-10);
       
    % updata theta
    tmp1 = 1/trace(H'*(D1-A1)*H); tmp2 = 1/trace(H'*(D2-A2)*H); tmp = tmp1+tmp2;
    theta1 = tmp1/tmp; theta2 = tmp2/tmp;
    
    if stop_rule == 2
        obj = compute_obj(X1,X2,W1,W2,H,D1,A1,D2,A2,gamma,theta1,theta2);
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

end

function obj = compute_obj(X1,X2,W1,W2,H,D1,A1,D2,A2,gamma,theta1,theta2)
  L1 = D1-A1; L2 = D2-A2;
  obj = norm(X1-W1*H','fro')^2+norm(X2-W2*H','fro')^2+gamma*(trace(theta1.^2*H'*L1*H+theta2.^2*H'*L2*H));
end
