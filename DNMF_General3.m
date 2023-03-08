function [cROIs, Cs, coherence, skew, sz] = DNMF_General3(V, options)
    % [cROIs, Cs, coherence, skew, sz] = DNMF_General3(V, options)
    
    %Extract options
    thr = options.thr;
    patchSize = options.patchSize;
    stride = options.stride;
    overlapThr = options.overlapThr;
    sizeRange = options.sizeRange;
    filtSize = options.filtSize;
    DETREND_FRAMES = options.DETREND_FRAMES;
    shapeThr = options.shapeThr;
    med_opt = options.med_opt;
    height = size(V,1);
    width = size(V,2);

    CORE_OR_RANDOM = 1;
    %% Define ROI Pieces & Activity Traces
    
    % Split video into overlapping spatial patches
    indicesA = 1:stride:(height-patchSize(1))+1;
    indicesB = 1:stride:(width-patchSize(2))+1;

    ROIs = cell(1, length(indicesA)*length(indicesB));
    Cs = cell(length(indicesA)*length(indicesB), 1);
    COHERE = cell(length(indicesA)*length(indicesB), 1);
    SKEW = cell(length(indicesA)*length(indicesB), 1);
    SIZE = cell(length(indicesA)*length(indicesB), 1);

    stamp = zeros(height, width);
    count = 1;
    % For each patch 
    for i_A = indicesA
        aa = i_A:(i_A+patchSize(1)-1);
        for i_B = indicesB
            bb = i_B:(i_B+patchSize(2)-1);
            fprintf('%d-%d\n',i_A,i_B);
            thisV = double(V(aa,bb,:));

            % Detrend this patch
            dv = j_detrend2b(thisV,DETREND_FRAMES,3,true,med_opt); 
    
            medThisV = nanmedian(dv.*zeroOut(dv>=0),3);
            dv(isnan(dv) | isinf(dv)) = 0;
            if(CORE_OR_RANDOM==1)
                % Find cores from original Y 
                A0 = findCores(dv, thr*medThisV, [1 1 filtSize],min(sizeRange)/2);
                A0 = reshape(A0,patchSize(1)*patchSize(2),[]);
                a0 = sparse(A0);
                roiAND = a0'*a0;
                roiOR = prod(patchSize) - (1-a0)'*(1-a0);
                roiJAC = roiAND./roiOR;
                [~,groups] = graphconncomp(roiJAC>overlapThr);
                temp = zeros(size(A0,1),max(groups));
                for i_group = 1:max(groups)
                    temp(:,i_group) = sum(A0(:,groups==i_group),2)>0;
                end
                A0 = temp;
            else
                % or Initialize Randomly
                A0 = 50;
            end
            
            defoptions = CNMFSetParms;
            defoptions.d1 = size(thisV,1);
            defoptions.d2 = size(thisV,2);
            defoptions.block_size = patchSize;
            defoptions.min_pixel = sizeRange(1);
            defoptions.thr_method = options.thr_method;
            defoptions.maxthr = options.maxthr;
            defoptions.quantileThr = options.quantileThr;
            defoptions.final_C = options.final_C;
            defoptions.eta = options.eta;
            defoptions.beta = options.beta;
            if(isempty(A0))
                continue;
            end
            % Pass through CNMF for A, C estimates
            [A,C,B,F] = sparse_NMF_initialization7(thisV,A0,defoptions);
            if(isempty(A))
                continue;
            end
            
            [coherence, skew, roiSize] = evaluateROIs(A, C, patchSize);
            
            % Calculate the correlation of the ROI with the video at the
            % frame of maximum activity, threshold based on that
            [~,i] = max(C,[],2);
            ppp = corrTrace(thisV,reshape(full(A),[patchSize(1) patchSize(2) size(A,2)]),3);
            w = arrayfun(@(x) ppp(x,i(x)),1:size(i));
            valid = w>=shapeThr;
            
            A = A(:,valid);
            C = C(valid,:);
            coherence = coherence(valid);
            skew = skew(valid);
            roiSize = roiSize(valid);
            
            % Plotting
            temp = zeros(height, width, size(A,2));
            temp2 = zeros(height, width, size(A,2));           
            AA = reshape(bsxfun(@times,full(A),max(C,[],2)'), patchSize(1), patchSize(2), size(A,2));
            AA2 = reshape(full(A), patchSize(1), patchSize(2), size(A,2));
            
            temp(aa,bb,:) = AA;        
            temp2(aa,bb,:) = AA2;
            ROIs{count} = sparse(reshape(temp2,height*width,[]));
            Cs{count} = full(C);
            COHERE{count} = coherence(:);
            SKEW{count} = skew(:);
            SIZE{count} = roiSize(:);
            count = count+1;
            if(~isempty(temp))
                stamp = max(stamp, max(temp,[],3));
                clf;
                subplot(1,2,1);
                imagesc(stamp);
                axis square;
                subplot(2,2,2);
                imagesc(max(thisV,[],3));
                axis square;
                subplot(2,2,4);
                imagescc(max(A*C,[],2));
                axis square;
                
                drawnow;
            end
        end
    end
    
    % Collect across patches
    fprintf('Collecting ROIs...');
    cROIs = cell2mat(ROIs);
    Cs = cell2mat(Cs);
    coherence = cell2mat(COHERE);
    skew = cell2mat(SKEW);
    sz = cell2mat(SIZE);
    fprintf('Done!\n');
end