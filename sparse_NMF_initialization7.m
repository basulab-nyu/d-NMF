function [Ain,Cin,bin,fin] = sparse_NMF_initialization7(Y,A0,options) %beta,eta,X0,err_thr,max_iter)

T = size(Y,ndims(Y));
if ~ismatrix(Y)
    Y = reshape(Y,numel(Y)/T,T);
end

defoptions = CNMFSetParms;
if nargin < 3 || isempty(options)
    options = defoptions; 
end
if ~isfield(options,'snmf_max_iter')
    options.snmf_max_iter = defoptions.snmf_max_iter; 
end
max_iter = options.snmf_max_iter;
if ~isfield(options,'err_thr') 
    options.err_thr = defoptions.err_thr; 
end
err_thr = defoptions.err_thr;
if ~isfield(options,'eta')
    options.eta = defoptions.eta; 
end
eta = options.eta*max(Y(:))^2;
if ~isfield(options,'beta')
    options.beta = defoptions.beta; 
end
beta = options.beta;
if ~isfield(options,'nb') 
    options.nb = defoptions.nb; 
end 
nb = options.nb;
if ~isfield(options,'rem_prct') || isempty(options.rem_prct)
    options.rem_prct = defoptions.rem_prct; 
end
    
options.conn_comp = false;
% options.maxthr = 0.25;
repeat = 1;
iter = 1;
obj_ = 1e-10;

v = ver;
flag_optim = any(strcmp('Optimization Toolbox', {v.Name})); % check if optimization toolbox is present
if flag_optim
    if verLessThan('optim','6.3')
        min_options = optimset('Algorithm','interior-point','GradObj','On','Display','Off');
    else
        min_options = optimoptions('fmincon','Algorithm','interior-point','GradObj','On','Display','Off');
    end
end

% remove median
medY = prctile(Y,options.rem_prct,2);
Y = bsxfun(@minus, Y, medY);

if(numel(A0)==1)
    K = A0;
    C = rand(A0,T);
    A = max((Y*C')*pinv(C*C'+beta*ones(K,K)),0);
else    
    extra = 2;
%     A0 = [A0 rand(size(A0,1),extra)<0.1];
%     A1 = [A0 sum(A0,2)==0 ones(size(A0,1),1)];
    if(any(sum(A0)==size(A0,1)))
        A1 = [A0 sum(A0,2)==0];
    else
        A1 = [A0 sum(A0,2)==0 ones(size(A0,1),1)];
    end
    valid = sum(A1)>0;
    A1 = A1(:,valid);
%     A1 = A0;
    A = bsxfun(@rdivide, A1, sum(A1));
%     K0 = size(A0,2);
    K = size(A,2);
    C = max((A'*A + eta*eye(K))\(A'*Y),0);
%     C = max((A'*A)\(A'*Y),0);   
%     C = [C; max(C(:))*rand(extra,size(C,2))/2];
end

% A = A0;
flagBreak = false;
while (iter <= max_iter) && repeat
    A = (A + max((Y*C')*pinv(C*C'+beta*ones(K,K)),0))/2;
    C = (C + max((A'*A + eta*eye(K))\(A'*Y),0))/2;    
    
%     A = max((A + (Y*C')*pinv(C*C'+beta*ones(K,K)))/2,0);
%     C = max((C + (A'*A + eta*eye(K))\(A'*Y))/2,0);
    
    ff = find(sum(C,2)==0);
    if ~isempty(ff)
        A(:,ff) = [];
        C(ff,:) = [];
        K = K - length(ff);
    end
    
    iter = iter + 1;
    if mod(iter,10) == 0
%         A = threshold_components2(A,options);      
%         [A,C] = component_split(A,C,options.block_size,20,4000);
%         K = size(A,2);
        if(isempty(A))
            flagBreak = true;
            break;
        end
        nC = double(sum(C.^2,2));
        AA = A'*A;
        mine = @(e) min_e(e,double(beta),double(eta),nC,A,AA);
        
        e = double(eta*nC./full(beta*A'*sum(A,2))).^(1/4);
%         if flag_optim && size(C,1)<100
%             e = fmincon(mine,max(e,1e-4),[],[],[],[],1e-4*ones(K,1),[],[],min_options);            
%         end
        C = diag(e)\C;
        A = A*diag(e);
%         fprintf('%i out of maximum %i iterations done \n',iter,max_iter);
        
        obj = norm(Y - A*C,'fro')^2 + eta*norm(C,'fro')^2 + beta*norm(sum(A,2))^2;
        repeat = abs(obj - obj_) > err_thr*obj_;
        obj_ = obj;
    end
    
end
% if(~flagBreak)

%     options.thr_method = 'quant';
%     options.quantileThr = 0.9;
    options.conn_comp = false;
    A1 = A;
    A = threshold_components3(A,options);
    [A,C] = component_split(A,C,options.block_size,options.min_pixel,Inf);
    K = size(A,2);
    if(K==0)
        flagBreak = true;
        C = [];
    else
        if(options.final_C)
            C = max((A'*A + eta*eye(K))\(A'*Y),0);
        end
    end
% end
fprintf('Algorithm converged after %i iterations. \n',iter-1);
% A = max((Y*C')*pinv(C*C'+beta*ones(K,K)),0);
% C = max((A'*A + eta*eye(K))\(A'*Y),0);    
% 
% ff = find(sum(C,2)==0);
% if ~isempty(ff)
%     A(:,ff) = [];
%     C(ff,:) = [];
%     K = K - length(ff);
% end
% options.conn_comp = false;
% options.medw = [1 1];
% A = threshold_components(A,options);       
% nC = double(sum(C.^2,2));
% AA = A'*A;
% mine = @(e) min_e(e,double(beta),double(eta),nC,A,AA);
% 
% e = double(eta*nC./full(beta*A'*sum(A,2))).^(1/4);
% if flag_optim
%     e = fmincon(mine,max(e,1e-4),[],[],[],[],1e-4*ones(K,1),[],[],min_options);            
% end
% C = diag(e)\C;
% A = A*diag(e);

        
        
        
        
        
[Ain,Cin] = order_components(A,C);
if(flagBreak)
    bin = [];
    fin = [];
else
    [bin,fin] = nnmf(max(Y - Ain*Cin + repmat(medY,1,T),0),nb);
end
Ain = sparse(double(Ain));
bin = double(bin);
    function [f,grad] = min_e(e,beta,eta,nC,A,AA)
        f = eta*norm(sqrt(nC)./e)^2 + beta*norm(A*e)^2;
        grad = -2*eta*nC./(e.^3) + 2*beta*AA*e;
    end


    function [A_or,C_or] = order_components(A,C)
        nA = sqrt(sum(A.^2));
        nr = length(nA);
        A = bsxfun(@times,A,1./nA(:)'); %A/spdiags(nA(:),0,nr,nr);
        mA = max(A);
        C = bsxfun(@times,C,nA(:)); %spdiags(nA(:),0,nr,nr)*C;
        nC2 = sqrt(sum(C.^2,2));
        mC = max(C,[],2);
        [~,srt] = sort(mC.*mA'./nC2,'descend');
        A_or = A(:,srt);
        C_or = C(srt,:);
    end

end