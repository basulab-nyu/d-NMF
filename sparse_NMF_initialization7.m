function [Ain,Cin,bin,fin] = sparse_NMF_initialization7(Y,A0,options) %beta,eta,X0,err_thr,max_iter)
    % [Ain,Cin,bin,fin] = sparse_NMF_initialization7(Y,A0,options) %beta,eta,X0,err_thr,max_iter)
    T = size(Y,ndims(Y));
    if ~ismatrix(Y)
        Y = reshape(Y,numel(Y)/T,T);
    end

    % Extract options
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
    repeat = 1;
    iter = 1;
    obj_ = 1e-10;

    
    % remove median
    medY = prctile(Y,options.rem_prct,2);
    Y = bsxfun(@minus, Y, medY);

    
    if(numel(A0)==1)    % Random initialization
        K = A0;
        C = rand(A0,T);
        A = max((Y*C')*pinv(C*C'+beta*ones(K,K)),0);
    else                % Seeded initialization    
        if(any(sum(A0)==size(A0,1)))
            A1 = [A0 sum(A0,2)==0];
        else
            A1 = [A0 sum(A0,2)==0 ones(size(A0,1),1)];
        end
        valid = sum(A1)>0;
        A1 = A1(:,valid);
        A = bsxfun(@rdivide, A1, sum(A1));
        K = size(A,2);
        C = max((A'*A + eta*eye(K))\(A'*Y),0);
    end

    flagBreak = false;
    while (iter <= max_iter) && repeat
        A = (A + max((Y*C')*pinv(C*C'+beta*ones(K,K)),0))/2;
        C = (C + max((A'*A + eta*eye(K))\(A'*Y),0))/2;    

        ff = find(sum(C,2)==0);
        if ~isempty(ff)
            A(:,ff) = [];
            C(ff,:) = [];
            K = K - length(ff);
        end

        iter = iter + 1;
        if mod(iter,10) == 0
            if(isempty(A))
                flagBreak = true;
                break;
            end
            nC = double(sum(C.^2,2));

            e = double(eta*nC./full(beta*A'*sum(A,2))).^(1/4);
            C = diag(e)\C;
            A = A*diag(e);

            obj = norm(Y - A*C,'fro')^2 + eta*norm(C,'fro')^2 + beta*norm(sum(A,2))^2;
            repeat = abs(obj - obj_) > err_thr*obj_;
            obj_ = obj;
        end

    end
    
    options.conn_comp = false;
    
    % Threshold and split ROIs
    A = threshold_components3(A,options);
    [A,C] = component_split(A,C,options.block_size,options.min_pixel,Inf);
    K = size(A,2);
    if(K==0)
        flagBreak = true;
        C = [];
    else
        if(options.final_C)
            C = max((A'*A + eta*eye(K))\(A'*Y),0);      %Re-compute C one more time
        end
    end

    fprintf('Algorithm converged after %i iterations. \n',iter-1);
        
    [Ain,Cin] = order_components(A,C);
    if(flagBreak)
        bin = [];
        fin = [];
    else
        [bin,fin] = nnmf(max(Y - Ain*Cin + repmat(medY,1,T),0),nb);
    end
    Ain = sparse(double(Ain));
    bin = double(bin);
            

    
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