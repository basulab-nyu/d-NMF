function [cROIs, Cs, cROIs_BG, Cs_BG, coherence, skew, sz, patch_ID, A0] = DNMF_General(paths_to_videos, options)
    % [cROIs, Cs, cROIs_BG, Cs_BG, coherence, skew, sz, patch_ID, A0] = DNMF_General(paths_to_videos, options)
    %
    % Jason Moore, 2025
    
    
    %Extract options
    maxVal = options.maxVal;
    spatial_smooth = options.spatial_smooth;
    temporal_filter_corner = options.temporal_filter_corner;
    fs = options.fs;
    SNR_thr = options.SNR_thr;
    thr = options.thr;
    idealPatchSize = options.patchSize;
    stride = options.stride;
    overlapThr = options.overlapThr;
    sizeRange = options.sizeRange;
    filtSize = options.filtSize;
    DETREND_FRAMES = options.DETREND_FRAMES;
    shapeThr = options.shapeThr;
    med_opt = options.med_opt;
    height = options.szA;
    width = options.szA;

    if(~isfield(options,'loadInPatches'))
        loadInPatches = false;
    else
        loadInPatches = options.loadInPatches;
    end
    
    if(~isfield(options,'DOWNSAMPLE_FACTOR'))
        DOWNSAMPLE_FACTOR = 20;
    else
        DOWNSAMPLE_FACTOR = options.DOWNSAMPLE_FACTOR;
    end

    CORE_OR_RANDOM = 1;
    
    %%
    if(temporal_filter_corner>0)
        try
            [b_filt, a_filt] = butter(2,2*temporal_filter_corner/fs,'low');
            filter_opt = true;
        catch e
            fprintf('\nCorner frequency %.02f not compatible with sampling rate %.1f.\nSkipping temporal filter.\n',temporal_filter_corner,fs);
            b_filt = NaN;
            a_filt = NaN;
            filter_opt = false;
        end
    else
        b_filt = NaN;
        a_filt = NaN;
        filter_opt = false;
    end
    
    
    %% Define ROI Pieces & Activity Traces
    
    % Split video into overlapping spatial patches
    
    indicesA = 1+(idealPatchSize(1):stride:height+stride-1)-idealPatchSize(1);
    indicesB = 1+(idealPatchSize(2):stride:width+stride-1)-idealPatchSize(2);
    
    ROIs = cell(1, length(indicesA)*length(indicesB));
    Cs = cell(length(indicesA)*length(indicesB), 1);    
    COHERE = cell(length(indicesA)*length(indicesB), 1);
    SKEW = cell(length(indicesA)*length(indicesB), 1);
    SIZE = cell(length(indicesA)*length(indicesB), 1);
    PATCH_ID = cell(length(indicesA)*length(indicesB), 1);
    
    stamp = zeros(height, width);
    
    if(~loadInPatches)
        V = cell(1,1,length(paths_to_videos));
        for i_file = 1:length(paths_to_videos)
            temp = bigread2(paths_to_videos{i_file});
            temp(temp>maxVal) = 0;
            V{i_file} = temp;
        end
        V = cell2mat(V);
        Cs_BG = zeros(length(indicesA)*length(indicesB), size(V,3));
    end
    
    cROIs_BG = zeros(height*width, length(indicesA)*length(indicesB));
    

    % For each patch 
    for i_A = indicesA
        aa = i_A:min([(i_A+idealPatchSize(1)-1) height]);
        for i_B = indicesB
            bb = i_B:min([(i_B+idealPatchSize(2)-1) width]);
            count = (find(indicesA==i_A)-1)*length(indicesB)+find(indicesB==i_B);
            if(loadInPatches)
                V = cell(1,1,length(paths_to_videos));
                for i_file = 1:length(paths_to_videos)
                    V{i_file} = bigread2_partial(paths_to_videos{i_file}, aa, bb);
                end
                V = cell2mat(V);
                if(count==1)
                    Cs_BG = zeros(length(indicesA)*length(indicesB), size(V,3));
                end
            end
            
            fprintf('%d-%d\n',i_A,i_B);
            if(loadInPatches)
                thisV0 = double(V);
            else
                thisV0 = double(V(aa, bb, :));
            end          
            patchSize = [size(thisV0,1) size(thisV0,2)];
            
            if(filter_opt)
                thisV02 = reshape(thisV0,prod(patchSize),[])';
                f_thisV02 = filtfilt(b_filt,a_filt,thisV02);
                thisV0 = reshape(f_thisV02',[patchSize(1) patchSize(2) size(thisV0,3)]);
            else
%                 thisV0 = thisV0;
            end
            
            if(spatial_smooth>0)
                thisV0 = gaussian_smooth_2d_slices(thisV0,spatial_smooth);
            else
%                 thisV0 = thisV0;
            end            
            thisV02 = reshape(thisV0,prod(patchSize),[])';     
            
            % For initialization, we downsample by a factor of DOWNSAMPLE_FACTOR
            thisV02 = [thisV02; NaN(ceil(size(thisV02,1)/DOWNSAMPLE_FACTOR)*DOWNSAMPLE_FACTOR-size(thisV02,1),size(thisV02,2))];
            thisV02b = reshape(thisV02,DOWNSAMPLE_FACTOR,[],size(thisV02,2));
            thisV = reshape(squeeze(nanmean(thisV02b,1))',patchSize(1),patchSize(2),[]);
            
            % Detrend this patch
            [dv,baseline] = j_detrend2b(thisV,ceil(DETREND_FRAMES/DOWNSAMPLE_FACTOR),3,true,med_opt); 
    
            medThisV = nanmedian(dv.*zeroOut(dv>=0),3);
            dv(isnan(dv) | isinf(dv)) = 0;
            if(CORE_OR_RANDOM==1)
                % Find cores from original Y 
                A0 = findCores(dv, thr*medThisV, [1 1 filtSize],min(sizeRange));
                if(isempty(A0))
                    continue;
                end
                A0 = reshape(A0,[],size(A0,3));
                a0 = sparse(A0);
                roiAND = a0'*a0;
                roiOR = prod(patchSize) - (1-a0)'*(1-a0);
                roiJAC = roiAND./roiOR;
                groups = conncomp(graph(roiJAC>overlapThr));
                temp = zeros(size(A0,1),max(groups));
                for i_group = 1:max(groups)
                    temp(:,i_group) = sum(A0(:,groups==i_group),2)>0;
                end
                A0 = temp;
                
                cThisV = reshape(thisV0,[],size(thisV0,3));
                C0 = A0'*cThisV;
                C0 = j_detrend2b(C0,DETREND_FRAMES,2,false,true);
                SNR = quantile(C0,0.999,2)./mad(C0,[],2);
                valid = SNR>SNR_thr;
                A0 = A0(:,valid);
                C0 = C0(valid,:);
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
            [A,C,B,F] = NMF_adapt3(thisV0,A0,defoptions);
            if(isempty(A))
                continue;
            end
            
            [coherence, skew, roiSize] = evaluateROIs(A, C, patchSize);
            
            % Calculate the correlation of the ROI with the video at the
            % frame of maximum activity, threshold based on that
            [~,i] = max(C,[],2);
            ppp = corrTrace(thisV0,reshape(full(A),[patchSize(1) patchSize(2) size(A,2)]),3);
            w = arrayfun(@(x) ppp(x,i(x)),1:length(i));
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

            temp2 = zeros(height, width);           
            temp2(aa,bb) = reshape(B, patchSize(1), patchSize(2));
            cROIs_BG(:,count) = temp2(:);
            Cs{count} = full(C);
            Cs_BG(count,:) = F;
            COHERE{count} = coherence(:);
            SKEW{count} = skew(:);
            SIZE{count} = roiSize(:);
            PATCH_ID{count} = count*ones(size(roiSize(:)));

            if(~isempty(temp))
                stamp = max(stamp, max(temp,[],3));
                clf;
                subplot(1,2,1);
                imagesc(stamp);
                hold on;
                plot([min(bb) max(bb)], min(aa)*[1 1], 'r', 'Linewidth', 2);
                plot([min(bb) max(bb)], max(aa)*[1 1], 'r', 'Linewidth', 2);
                plot(min(bb)*[1 1], [min(aa) max(aa)], 'r', 'Linewidth', 2);
                plot(max(bb)*[1 1], [min(aa) max(aa)], 'r', 'Linewidth', 2);                
                axis square;
                xlabel('X (px)');
                ylabel('Y (px)');
                
                subplot(2,2,2);
                imagesc(min(bb)+(0:size(thisV,2)-1), min(aa)+(0:size(thisV, 1)-1), max(thisV, [], 3));
                axis square;
                xlabel('X (px)');
                ylabel('Y (px)');
                
                subplot(2,2,4);
                imagesc(min(bb)+(0:size(thisV,2)-1), min(aa)+(0:size(thisV, 1)-1), reshape(max(A*C, [], 2), [patchSize(1) patchSize(2)]));
                axis square;
                xlabel('X (px)');
                ylabel('Y (px)');
                
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
    patch_ID = cell2mat(PATCH_ID);
    fprintf('Done!\n');
end