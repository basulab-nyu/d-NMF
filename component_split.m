function [A2, C2] = component_split(A, C, shape, minSize, maxSize)

    if(nargin<5)
        maxSize = Inf;
    end
    if(nargin<4)
        minSize = 1;
    end
    

    nComponents = size(A,2);

    
    A1 = cell(1,nComponents);
    C1 = cell(nComponents,1);
    for i_component = 1:nComponents
        this = reshape(A(:,i_component), shape);
        mask = this>0;
        CC = bwconncomp(full(mask));
        sz = cellfun(@(x) length(x), CC.PixelIdxList);
        valid = isbetween(sz, minSize, maxSize);
                
        if(sum(valid)>0)

            pixelLists = CC.PixelIdxList(valid);
            temp = zeros(size(this));
            outA2 = zeros(size(A,1),length(pixelLists));
            for i_sub_component = 1:length(pixelLists)
                temp(pixelLists{i_sub_component}) = i_sub_component;
                outA2(:,i_sub_component) = reshape(this.*(temp==i_sub_component),[],1);
            end
            A1{i_component} = outA2;
            C1{i_component} = repmat(C(i_component,:),length(pixelLists),1);
        end
    end
    A2 = cell2mat(A1);
    C2 = cell2mat(C1);
    
    
end