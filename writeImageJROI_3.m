% Writes all features of the ImageJ (IJ) binary roi file format (*.roi) as specified by ij.io.RoiDecoder.
%    matlabpathcoors - an ordered list of coordinate pairs that describe a single ROI boundary.  Coordinates assumed to be (y,x) order (where y is vertical image axis and x horizontal image axis) and follow the Matlab “origin @ (1,1)” convention.
%    method - The IJ 'roi_type' that will be used.
        % Possiblities are:
        % 1 - Use 'freehandline' type.
        % 2 - Use 'freehandselection' type (unclosed shapes are automatically closed by ImageJ).
        % 3 - Coordinate lists where first and last coordinate are within 1 pixel's distance are automatically set as 'freehandselection' type, else set to 'freehandline' type
%    ch,z,t - integers representing the postion of the roi within a larger hyperstack. Set ch,z,t = 0 when image is a single plane.
%    folder – A file system path in char array format pointing to where the ImageJ binary roi file will be written.  Must end in '\'.
        % An ImageJ-style roi filename will be automatically generated and append to this path.

function writeImageJROI_3(matlabpathcoors, method, z, roi_name, folder)

    ch = 0;
%     z = 0;
    t = z;
    
    % -- A table array template that holds default parameters, describes data types, and describes the correspondence btwn matlab and java data
    % ROI-specific data is added to this table.
    % See Matlab_writeIJROI_Function_Description.docx for a complete description of ROI_FILE_DATA.mat
    load('ROI_FILE_DATA');

    % -- Processes a Matlab coordinate list into the IJ ROI list format
    % Eg. finds bounding box, converts to relative coordinate system, puts coordinates in proper order  
    % Returns derivative ij roi information such as number of coordinates, header2_offset, etc.
    [roi_type, box_top, box_left, box_bottom, box_right, num_roi_coors, header2_offset, roi_coor_list] = getroidata(matlabpathcoors, method);

    % -- Assigns above roi data to proper cells in ROI_FILE_DATA
    ROI_FILE_DATA.integer_data{'roi_type'} = roi_type;
    ROI_FILE_DATA.integer_data{'box_top'} = box_top;
    ROI_FILE_DATA.integer_data{'box_left'} = box_left;
    ROI_FILE_DATA.integer_data{'box_bottom'} = box_bottom;
    ROI_FILE_DATA.integer_data{'box_right'} = box_right;
    ROI_FILE_DATA.integer_data{'num_roi_coors'} = num_roi_coors;
    ROI_FILE_DATA.integer_data{'header2_offset'} = header2_offset;
    ROI_FILE_DATA.integer_data{'roi_coor_list'} = roi_coor_list;

    % -- Accessory information for rois positioned within a hyperstack
    ROI_FILE_DATA.integer_data{'ch_position'} = ch;
    ROI_FILE_DATA.integer_data{'z_position'} = z;
    ROI_FILE_DATA.integer_data{'t_position'} = t;

    %  -- Creates an ROI name according to the ij default style
    slice = 1; % Slice = 1 works for also for single slice images. IJ will just ignore the slice value.
    cntry = round((box_bottom - box_top)/2)+box_top; 
    cntrx = round((box_right - box_left)/2)+box_left;
%     roi_name = getroiname(slice, cntry, cntrx); % roiname is an ij standard format roiname string represented as a matlab character array
    % derivative information
    name_length = length(roi_name); % in text symbols
    name_offset = header2_offset + 64; % header2 is 64 bytes long

    % Creates full file path
    filepath = [folder,roi_name,'.roi']; % appends char arrays

    % -- Converts roi_name into a 16-bit character encoding format as need to write the binary
    % This format is not directly supported in matlab so this is a kluge
    long_name_length = 2*name_length; % x2 because ij roi format uses 2 bytes per symbol
    long_roi_name = char(zeros(1,long_name_length)); % initializes to all null characters (ascii corresponding to 0)
    for i = 1:length(roi_name)
        % writes a roi_name character symbol (non-null) to every other index (resulting in 16 bits per char).
        long_roi_name(i*2) = roi_name(i);
    end

    % -- Assigns roi name related data
    ROI_FILE_DATA.integer_data{'roi_name'} = long_roi_name; % two bytes per char with first byte of each pair always null
    ROI_FILE_DATA.integer_data{'name_offset'} = name_offset;
    ROI_FILE_DATA.integer_data{'name_length'} = name_length; % number of text symbols

    % -- Writes the information in ROI_FILE_DATA as a binary according to IJ standards
    permission = 'w'; % for write
    machinefmt = 'b'; % big-endian per ij
    encodingIn = 'US-ASCII'; % character encoding. Does not matter here because all data is integers written as byte(s).
    fileID = fopen(filepath,permission,machinefmt,encodingIn);
    % Writes data sequentially to the fileID
    writeroidata(ROI_FILE_DATA,fileID,machinefmt);

end % EO function



% -----  Supporting Functions -------

function [roi_type, box_top, box_left, box_bottom, box_right, num_roi_coors, header2_offset, roi_coor_list] = getroidata(matlabpathcoors, method)
    
    % Converts to ij origin style (0,0) rather than (1,1)
    ijpathcoors = matlabpathcoors - 1;

    % Method sets how roi_type is determined.  Methods are:
    % 1 - Use freehandline for all
    % 2 - Use freehandselection for all (are unclosed shapes automatically closed by IJ???)
    % 3 - Closed shapes are set as freehandselection, open shapes as freehandline
    if method == 1
        roi_type = 4; % freeline, most common
    end
    if method == 2
        roi_type = 7; % freehand selection
    end
    if method == 3
        % Tests for closed shape vs open shape (line). If first and last pix are adjacent, it is a closed shape
        end2enddist = pdist2(ijpathcoors(1,:),ijpathcoors(end,:),'euclidean');
        if end2enddist < 1.5 % accounts for diagonal adjacency
            roi_type = 7;
        else
            roi_type = 4;
        end
    end

    if method == 4
        roi_type = 0;
    end
    % Finds edges of bounding box.
    ypathcoors = ijpathcoors(:,1);
    xpathcoors = ijpathcoors(:,2);

    box_top = min(ypathcoors,[],1);
    box_bottom = max(ypathcoors,[],1);
    box_left = min(xpathcoors,[],1);
    box_right = max(xpathcoors,[],1);
    % fixes the problem that min and max return ALL values
    box_top = box_top(1);
    box_bottom = box_bottom(1);
    box_left = box_left(1);
    box_right = box_right(1);

    % Finds number of coordinates (pairs)
    num_roi_coors = size(ijpathcoors,1);

    % Finds number of bytes to start of header2. Bytes in first header plus bytes in coor list.
    header2_offset = 64+(num_roi_coors*2*2); % 2 values per coordinate x 2 bytes per value (signed short)

    % Creates IJ style roi coordinate list
    dycoors = ypathcoors - box_top;
    dxcoors = xpathcoors - box_left;
    roi_coor_list = [dxcoors; dycoors]; % appends column vectors
    roi_coor_list = transpose(roi_coor_list); % row vector is assumed when writting below...
end


function [roi_name] = getroiname(slice, cntry, cntrx)
    % Per the IJ ROI Manager standard, the ROI name is three blocks of four numbers ZZZZ-YYYY-XXX where zzzz is the slice number, yyyy is the center coordinate of the ROI bounding box in Y and xxxx is the center coordinate of the ROI bounding box in X. 
    % Inputs assumed to be appropriate positive integers

    zval = num2str(slice,'%04d'); % formatting operator produces 4 digit values
    yval = num2str(cntry,'%04d'); 
    xval = num2str(cntrx,'%04d'); 

    namestring = strcat(zval,'-',yval,'-',xval);
    roi_name = char(namestring); % need a char vector because each letter will be written separately.

end

function writeroidata(ROI_FILE_DATA,fileID,machinefmt)
    
    skip = 0; % don't skip bytes

    % Loops over ROI file data and writes it as binary
    % Data is inherently written sequentially to the file as long as the fileID is not closed
    for i = 1:size(ROI_FILE_DATA,1)
        
        swch = 0; % a logic control constant        

        % Special handling of rows where data itself is an array and so requires further looping to write
        if strcmp(ROI_FILE_DATA.matlab_var{i},'Iout')
            Iout_data = ROI_FILE_DATA.integer_data{'Iout'};
            for j = 1:size(Iout_data,2)
                fwrite(fileID,Iout_data(1,j),ROI_FILE_DATA.matlab_type{i},skip,machinefmt);
            end
            swch = 1;
        end
        if strcmp(ROI_FILE_DATA.matlab_var{i},'size')
            size_data = ROI_FILE_DATA.integer_data{'size'};
            for j = 1:size(size_data,2)
                fwrite(fileID,size_data(1,j),ROI_FILE_DATA.matlab_type{i},skip,machinefmt);
            end
            swch = 1;
        end
        if strcmp(ROI_FILE_DATA.matlab_var{i},'roi_coor_list')
            roi_coor_list_data = ROI_FILE_DATA.integer_data{'roi_coor_list'};
            for j = 1:size(roi_coor_list_data,2)
                fwrite(fileID,roi_coor_list_data(1,j),ROI_FILE_DATA.matlab_type{i},skip,machinefmt);
            end
            swch = 1;
        end
        if strcmp(ROI_FILE_DATA.matlab_var{i},'gap2')
            gap2_data = ROI_FILE_DATA.integer_data{'gap2'};
            for j = 1:size(gap2_data,2)
                fwrite(fileID,gap2_data(1,j),ROI_FILE_DATA.matlab_type{i},skip,machinefmt);
            end
            swch = 1;
        end
        if strcmp(ROI_FILE_DATA.matlab_var{i},'roi_name')
            roi_name_data = ROI_FILE_DATA.integer_data{'roi_name'};
            for j = 1:size(roi_name_data,2)
                fwrite(fileID,roi_name_data(1,j),ROI_FILE_DATA.matlab_type{i},skip,machinefmt);
            end
            swch = 1;
        end
        
        % Writting of rows where data is a scalar
        % swch ensures this only occurs if data was not already written above...
            % i.e. data is only written once for each row 
        if swch == 0
            fwrite(fileID,ROI_FILE_DATA.integer_data{i},ROI_FILE_DATA.matlab_type{i},skip,machinefmt);
        end

    end

    fclose(fileID);

end
