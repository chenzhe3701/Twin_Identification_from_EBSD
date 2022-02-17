
%% setup
clear; clc; close all;
[f_setting,p_setting] = uigetfile('','select setting file (.m format)');

% run the .m file as setting file to load variables
run(fullfile(p_setting,f_setting));
% we have these variables
assert(exist('working_dir','var')==1);
assert(exist('sample_name','var')==1);
assert(exist('iE_max','var')==1);
assert(exist('grain_pair','var')==1);
assert(exist('ID_list','var')==1);
assert(exist('ID_merge_list','var')==1);

save_dir = fullfile(working_dir, 'analysis');
mkdir(save_dir);
%% Step-1: load .txt grain files, align_euler_to_sample and save as .mat

% (1.1) regular grain file
for iE = 0:iE_max

    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,[sample_name,' grain_file_type_1 iE=',num2str(iE),'.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,[sample_name,' grain_file_type_2 iE=',num2str(iE),'.txt']));
    % find column
    column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1, ...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
    column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2, ...
        {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

    gID = EBSD_data_2(:,column_index_2(1));
    gPhi1 = EBSD_data_2(:,column_index_2(2));
    gPhi = EBSD_data_2(:,column_index_2(3));
    gPhi2 = EBSD_data_2(:,column_index_2(4));
    gCenterX = EBSD_data_2(:,column_index_2(5));
    gCenterY = EBSD_data_2(:,column_index_2(6));
    gEdge = EBSD_data_2(:,column_index_2(10));

    gNNeighbors = EBSD_data_2(:,column_index_2(7));
    gDiameter = EBSD_data_2(:,column_index_2(8));
    gArea = EBSD_data_2(:,column_index_2(9));
    gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));

    x = EBSD_data_1(:,column_index_1(5));
    y = EBSD_data_1(:,column_index_1(6));
    unique_x = unique(x(:));
    ebsdStepSize = unique_x(2) - unique_x(1);
    mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
    nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

    ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    phi1 = reshape(EBSD_data_1(:,column_index_1(2)),mResize,nResize)';
    phi = reshape(EBSD_data_1(:,column_index_1(3)),mResize,nResize)';
    phi2 = reshape(EBSD_data_1(:,column_index_1(4)),mResize,nResize)';
    x = reshape(EBSD_data_1(:,column_index_1(5)),mResize,nResize)';
    y = reshape(EBSD_data_1(:,column_index_1(6)),mResize,nResize)';
    % change it to degrees, if necessary
    if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
        phi1 = phi1*180/pi();
        phi = phi*180/pi();
        phi2 = phi2* 180/pi();
    end

    [phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 'non-mtex', 90,180,0);
    [gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'non-mtex', 90,180,0);

    euler_aligned_to_sample = 1;
    save(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

end

% (1.2) parent grain file. [[[For iE=0, use grain_file as parent_grain_file]]]
for iE = 0:iE_max
    [EBSD_data_1, EBSD_header_1] = grain_file_read(fullfile(working_dir,[sample_name,' parent_grain_file_type_1 iE=',num2str(iE),'.txt']));
    [EBSD_data_2, EBSD_header_2] = grain_file_read(fullfile(working_dir,[sample_name,' parent_grain_file_type_2 iE=',num2str(iE),'.txt']));
    % find column
    column_index_1 = find_variable_column_from_grain_file_header(EBSD_header_1, ...
        {'grain-ID','phi1-r','phi-r','phi2-r','x-um','y-um','edge'});
    column_index_2 = find_variable_column_from_grain_file_header(EBSD_header_2, ...
        {'grainId','phi1-d','phi-d','phi2-d','x-um','y-um','n-neighbor+id','grain-dia-um','area-umum','edge'});

    gID = EBSD_data_2(:,column_index_2(1));
    gPhi1 = EBSD_data_2(:,column_index_2(2));
    gPhi = EBSD_data_2(:,column_index_2(3));
    gPhi2 = EBSD_data_2(:,column_index_2(4));
    gCenterX = EBSD_data_2(:,column_index_2(5));
    gCenterY = EBSD_data_2(:,column_index_2(6));
    gEdge = EBSD_data_2(:,column_index_2(10));

    gNNeighbors = EBSD_data_2(:,column_index_2(7));
    gDiameter = EBSD_data_2(:,column_index_2(8));
    gArea = EBSD_data_2(:,column_index_2(9));
    gNeighbors = EBSD_data_2(:,(column_index_2(7)+1):(size(EBSD_data_2,2)));

    x = EBSD_data_1(:,column_index_1(5));
    y = EBSD_data_1(:,column_index_1(6));
    unique_x = unique(x(:));
    ebsdStepSize = unique_x(2) - unique_x(1);
    mResize = (max(x(:)) - min(x(:)))/ebsdStepSize + 1;
    nResize = (max(y(:)) - min(y(:)))/ebsdStepSize + 1;

    ID = reshape(EBSD_data_1(:,column_index_1(1)),mResize,nResize)';
    phi1 = reshape(EBSD_data_1(:,column_index_1(2)),mResize,nResize)';
    phi = reshape(EBSD_data_1(:,column_index_1(3)),mResize,nResize)';
    phi2 = reshape(EBSD_data_1(:,column_index_1(4)),mResize,nResize)';
    x = reshape(EBSD_data_1(:,column_index_1(5)),mResize,nResize)';
    y = reshape(EBSD_data_1(:,column_index_1(6)),mResize,nResize)';
    % change it to degrees, if necessary
    if max(phi1(:))<7 && max(phi(:))<7 && max(phi2(:))<7
        phi1 = phi1*180/pi();
        phi = phi*180/pi();
        phi2 = phi2* 180/pi();
    end

    [phi1, phi, phi2] = align_euler_to_sample(phi1, phi, phi2, 'non-mtex', 90,180,0);
    [gPhi1, gPhi, gPhi2] = align_euler_to_sample(gPhi1, gPhi, gPhi2, 'non-mtex', 90,180,0);

    euler_aligned_to_sample = 1;
    save(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');
end
% For iE=0, use grain_file as parent_grain_file
% copyfile(fullfile(save_dir, [sample_name,'_grain_file_iE_0.mat']), fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']), 'f');

% save a copy after step 1
save_dir_1 = fullfile(save_dir, 'step-1');
mkdir(save_dir_1);
for iE = 0:iE_max
    copyfile(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir_1, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), 'f');
    copyfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir_1, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'f');
end

%% Step-2: Analyze parent grain file, affine transform to align maps of different load steps
% [optionally but suggested] plot ID map and manually select control IDs
% for alignment, and save to variable 'grain_pair', where each row contains
% the control IDs at each load step
run(fullfile(p_setting,f_setting));
for iE = 0:iE_max
    iB = iE + 1;
    d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    x = d.x;
    y = d.y;
    boundary = find_one_boundary_from_ID_matrix(ID);
    myplot(x,y,ID, boundary);
    % if have control IDs, label them on map
    try
       label_map_with_ID(x,y,ID,gcf,grain_pair(iB,:),'r',12,1); 
    end
    title(['iE=',num2str(iE)],'fontweight','normal');
    saveas(gcf, fullfile(working_dir,'analysis',['A1 original ID map iE=',num2str(iE),'.fig']));  
    print(fullfile(working_dir,'analysis',['A2 original ID map iE=',num2str(iE),'.tif']), '-dtiff');
    close;
    myplot(auto_grain(ID), boundary);
    title(['iE=',num2str(iE)],'fontweight','normal');
    print(fullfile(save_dir, ['A3 auto grain map iE=',num2str(iE),'.tif']),'-r300','-dtiff');
    close;
end

%% (2.1) transform and align each load step
close all;
for iE = 0:iE_max
    % iB=iE+1 for 0-based indexing, iB = iE for 1-based indexing
    iB = iE + 1;

    % reference, iE=0
    d = matfile(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    ID_0 = d.ID;
    x_0 = d.x;
    y_0 = d.y;
    % deformed iE
    d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    x = d.x;
    y = d.y;

    % load updated variable 'grain_pair'
    run(fullfile(p_setting,f_setting));
    try
        control_IDs = grain_pair([1,iB],:);
    catch
        control_IDs = [];
    end
    [tform, return_state] = align_ID_maps(ID_0, ID, control_IDs, x_0,y_0, x,y);
    tform_temp{iB} = tform; % store this in a cell

    % affine transform based on the determined 'tform'
    ID_iB_to_1 = interp_data(x,y,ID, x,y,tform.invert, 'interp', 'nearest');
    boundary_iB_to_1 = find_boundary_from_ID_matrix(ID_iB_to_1);
    ID_cell{iB} = ID_iB_to_1;
    gb_cell{iB} = boundary_iB_to_1;

    if iB==1
        boundary_all = boundary_iB_to_1;
    else
        boundary_all(boundary_iB_to_1>0) = iB * boundary_iB_to_1(boundary_iB_to_1>0);
    end
end

myplot(boundary_all);
title('boudary of transformed ID maps of all load steps','fontweight','normal');

%% (2.2) Select grains to be divided into more grains
close all;

% run the .m file as setting file to load variables
run(fullfile(p_setting,f_setting));
assert(exist('ID_list','var')==1);

for iE = []%0:iE_max
    iB = iE + 1;
    ID_iB_to_1 = ID_cell{iB};
    f = myplot(x,y,auto_grain(ID_iB_to_1)-6,boundary_all);
    title(['iE=',num2str(iE),', select grain to divide, double click to confirm, X to exist.'],'fontweight','normal');
    caxis([-5 -1]);
    % draw existing ones
    divide_list = ID_list{iB};
    for id = divide_list
        label_map_with_ID(x,y,ID_iB_to_1,gcf,id,'r',12,1);
    end

    while isvalid(f)
        try
            H = drawpoint;
            customWait(H);
            xy = round(H.Position);
            divide_list = [divide_list, ID_iB_to_1(xy(2),xy(1))];
        end
    end
    divide_list_cell{iB} = divide_list;

    disp(['Selected grains to divide, ID_list{',num2str(iB),'}:']);
    disp(divide_list_cell{iB});
end

%% (2.3) Select grains to be merged
% run the .m file as setting file to load variables
run(fullfile(p_setting,f_setting));
assert(exist('ID_merge_list','var')==1);

for iE = []%0:iE_max
    iB = iE + 1;
    ID_iB_to_1 = ID_cell{iB};
    f = myplot(x,y,auto_grain(ID_iB_to_1)-6,boundary_all);
    title(['iE=',num2str(iE),', select grain to merge, red=from, white=to. Double click to confirm, X to exist.'],'fontweight','normal');
    caxis([-5 -1]);
    % draw existing ones
    merge_list = ID_merge_list{iB};
    for ii = 1:size(merge_list,1)
        id_from = merge_list(ii,1);
        id_to = merge_list(ii,2);
        label_map_with_ID(x,y,ID_iB_to_1,gcf,id_from,'r',12,1);
        label_map_with_ID(x,y,ID_iB_to_1,gcf,id_to,'w',12,1);
    end

    while isvalid(f)
        try
            H = drawpoint('color','r');
            customWait(H);
            xy_from = round(H.Position);

            H = drawpoint('color','w');
            customWait(H);
            xy_to = round(H.Position);

            merge_list = [merge_list; ID_iB_to_1(xy_from(2),xy_from(1)),ID_iB_to_1(xy_to(2),xy_to(1))];
        end
    end
    merge_list_cell{iB} = merge_list;

    disp(['Selected grains to divide, ID_merge_list{',num2str(iB),'}:']);
    disp(merge_list_cell{iB});
end

%% (2.4) Merge, then draw mask to divide
% run the .m file as setting file to load variables
run(fullfile(p_setting,f_setting));
assert(exist('ID_list','var')==1);
assert(exist('ID_merge_list','var')==1);

for iE = 0:iE_max
    %% can run cell for each iE

    iB = iE + 1;

    % parent grain data for this load step
    d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID = d.ID;
    x = d.x;
    y = d.y;
    phi1 = d.phi1;
    phi = d.phi;
    phi2 = d.phi2;
    gID = d.gID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;
    euler_aligned_to_sample = d.euler_aligned_to_sample;

    % child grain data
    d = load(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
    ID_c = d.ID;
    gID_c = d.gID;

    % load the saved masks
    try
        load(fullfile(save_dir,[sample_name,'_mask_cell.mat']),'mask_cell','ID_updated_cell');
    catch
        mask_cell = cell(1,iE_max+1);
    end

    iN = 1;
    ID_updated = ID;
    next_gID = max(gID) + 1;

    for ii = 1:size(ID_merge_list{iB},1)
        ID_updated(ID_updated==ID_merge_list{iB}(ii,1)) = ID_merge_list{iB}(ii,2);
    end
    boundary = find_one_boundary_from_ID_matrix(ID_updated);
    boundary_new = boundary;
    % In the next part, we need to check if the added boundary segment just create one addtional grain.

    %% So, here, the function is to :
    % (1) select the grain, (2) calculate local misorientation map, (3) find a
    % line with the largest misorientation to divide the grain into two.

    N = length(ID_list{iB});
    while iN <= N
        try
            close(hf);
        end
        str = sprintf('iE=%d, iN=%d', iE, iN);
        disp(str);
        % example for debug, iE=1, grain 57
        ID_current = ID_list{iB}(iN);
        ind_local = ismember(ID_updated, ID_current); %ismember(ID, [ID_current,ID_neighbor]);
        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');

        ID_local = ID_updated(indR_min:indR_max, indC_min:indC_max);    % crop from ID_temp
        ID_c_local = ID_c(indR_min:indR_max, indC_min:indC_max);

        x_local = x(indR_min:indR_max, indC_min:indC_max);
        y_local = y(indR_min:indR_max, indC_min:indC_max);
        phi1_local = phi1(indR_min:indR_max, indC_min:indC_max);
        phi_local = phi(indR_min:indR_max, indC_min:indC_max);
        phi2_local = phi2(indR_min:indR_max, indC_min:indC_max);
        boundary_local = find_one_boundary_from_ID_matrix(ID_local);

        % for each pixel, find max misorientation within its four neighbors
        misorientation_max = calculate_max_local_misorientation_hcp(phi1_local, phi_local, phi2_local);

        % if not previously processed, process, save the mask
        % plot the max misorientation map, use mask to select the region where you
        % want to add the misorientation boundary as a new boundary
        if (length(mask_cell{iB})>=iN) && ~isempty(mask_cell{iB}{iN}) && length(mask_cell{iB}{iN}(:))==length(ID_local(:))
            mask = mask_cell{iB}{iN};
        else
            try
                map_t = -auto_grain(ID_c_local);
                map_t(boundary_local==1)=0;
                hf2 = myplotm(map_t, 'x',x_local, 'y',y_local);
                caxism([-5, 0]);
                cmap = colormap;
                cmap(end,:) = [1 1 1];
                colormap(cmap);
                % hf2 = myplotm(boundary_local,'x',x_local,'y', y_local);
                label_map_with_ID(x_local,y_local,ID_local, gcf, ID_current, 'r',12,1,0);
                title(['Draw mask to cover child grain boundaries to divide grain.',newline,'If child grain boundaries cannot be used, X to proceed.']);
                h = drawpolygon;
                customWait(h);
            catch
                hf2 = myplotm(map_t, 'x',x_local, 'y',y_local);
                caxism([-5, 0]);
                cmap = colormap;
                cmap(end,:) = [1 1 1];
                colormap(cmap);
                hf3 = myplotm(misorientation_max);
                caxism([3, 100]);
                set(gca,'xlim', get(gca,'xlim')+[-1 1]);    % need to do this for Matlab r2021b as it does not allow to drawpolygon outside of limits
                set(gca,'ylim', get(gca,'ylim')+[-1 1]);
                h = drawpolygon;
                customWait(h);
            end
            mask = h.createMask();
        end

        [ID_local_updated, need_to_redraw, need_to_use_misorientation_max] = ...
            divide_grain(ID_local, ID_c_local, mask, misorientation_max, next_gID);

        if need_to_redraw > 0
            mask_cell{iB}{iN} = [];     % delete bad mask
            disp(['iN=',num2str(iN),', need to redraw: ',num2str(need_to_redraw)]);
            warning(' More than one additional number of grains is created, you might want to do it again!');
        else
            ID_updated(indR_min:indR_max, indC_min:indC_max) = ID_local_updated;    % update ID_new
            ID_updated_cell{iB} = ID_updated;
            mask_cell{iB}{iN} = mask;
            save(fullfile(save_dir, [sample_name,'_mask_cell.mat']),'mask_cell','ID_updated_cell');

            disp(['iN=',num2str(iN),', need to redraw: ',num2str(need_to_redraw),', need to use misorientation_max: ',num2str(need_to_use_misorientation_max)])

            next_gID = next_gID + 1;    % update next_gID and iN
            iN = iN + 1;
        end
        try
            close(hf2);
            close(hf3);
        end
    end

    close all;
    gb = find_one_boundary_from_ID_matrix(ID_updated);
    myplot(ID_updated, gb);
    title('updated ID map');

    %% After adding more gb and get ID_updated, match ID_updated to ID_iE, update ID and grain info
    % ID_updated: modified from the added grain boundaries
    % ID: the target ID map at this iE

    % umPerDp = 1;    % micron per data point in EBSD data
    % data_in.umPerDp = umPerDp;
    data_in.symmetry = 'hcp';
    data_in.ID_target = ID;
    data_in.ID_temp = ID_updated;
    data_in.x = x;
    data_in.y = y;
    data_in.phi1 = phi1;
    data_in.phi = phi;
    data_in.phi2 = phi2;
    data_in.gID = gID;
    data_in.gPhi1 = gPhi1;
    data_in.gPhi = gPhi;
    data_in.gPhi2 = gPhi2;

    data_out = match_ID_map_and_summarize(data_in);

    ID = data_out.ID;   % the modified ID map
    gID = data_out.gID;
    gPhi1 = data_out.gPhi1;
    gPhi = data_out.gPhi;
    gPhi2 = data_out.gPhi2;
    gCenterX = data_out.gCenterX;
    gCenterY = data_out.gCenterY;
    gNNeighbors = data_out.gNNeighbors;
    gDiameter = data_out.gDiameter;
    gArea = data_out.gArea;
    gEdge = data_out.gEdge;
    gNeighbors = data_out.gNeighbors;

    save_dir_2 = fullfile(save_dir, 'step-2');
    mkdir(save_dir_2);

    save(fullfile(save_dir_2, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

    close all;

    disp(['finished iE=',num2str(iE)]);

    % try to copy reference immediately. But if need to redo, need to copy from
    % the saved file in step-1 folder.
    if iE == 0
        copyfile(fullfile(save_dir_2, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
            fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'f');
    end

end

%%
for iE = 0:iE_max
    copyfile(fullfile(save_dir_2, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'f');
end

%% [clean (child) grain file] For twin grains, if it contains more than one parent grain, then divide it
for iE = 0:iE_max
    %% do this iE
    iB = iE + 1;

    % parent grain file
    d = load(fullfile(save_dir_2, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    ID_parent = d.ID;

    % deformed iE
    d = load(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
    % read type-2 grain file and get average info for grains
    gID = d.gID;
    gPhi1 = d.gPhi1;
    gPhi = d.gPhi;
    gPhi2 = d.gPhi2;

    ID = d.ID;
    phi1 = d.phi1;
    phi = d.phi;
    phi2 = d.phi2;
    euler_aligned_to_sample = d.euler_aligned_to_sample;

    boundary = find_one_boundary_from_ID_matrix(ID);

    iN = 1;
    ID_updated = ID;     % up date child grain IDs
    next_gID = max(gID) + 1;

    boundary = find_one_boundary_from_ID_matrix(ID_updated);
    boundary_new = boundary;

    for ii = 1:length(gID)
        ID_current = gID(ii);
        ind_local = ismember(ID_updated, ID_current);

        indC_min = find(sum(ind_local, 1), 1, 'first');
        indC_max = find(sum(ind_local, 1), 1, 'last');
        indR_min = find(sum(ind_local, 2), 1, 'first');
        indR_max = find(sum(ind_local, 2), 1, 'last');

        ID_p_overlap = unique(ID_parent(ind_local));

        % if overlap with multiple parent grain (2021-07-17 note: some of these might be created due to the manual further division, as the boundary is difficult to control)
        if length(ID_p_overlap)>1
            disp(['===>',num2str(ii)]);
            vols = arrayfun(@(x) sum(ID_parent(ind_local)==x), ID_p_overlap); % overlap area with each of the overlapped parent ID
            [~, index_in_sorted] = sort(vols, 'descend');
            ID_p_overlap_sorted(index_in_sorted) = ID_p_overlap; % sort the overlapped parent ID, the one with the max overlap is sorted 1st

            for jj=2:length(ID_p_overlap_sorted)
                id_p = ID_p_overlap_sorted(jj);
                inds = ismember(ID_updated, ID_current) & ismember(ID_parent, id_p);

                % modify the overlapped part of this twin grain to have a new ID
                ID_updated(inds) = next_gID;
                % Note: 2022-01-04: an improved method may be required to prevent forming new child grains with disconnected pixels
                % Method: for inds area, first do one-pass-label, then assign each 'grain' with a new 'next_gID'
                % But still, these [new isolated small] grains should be removed from twin identification analysis

                % increment the next_gID
                next_gID = next_gID + 1;
            end

            % update the grain boundary [for illustration]
            ID_updated_local = ID_updated(indR_min:indR_max, indC_min:indC_max);
            boundary_local_new = find_one_boundary_from_ID_matrix(ID_updated_local);
            boundary_new(indR_min:indR_max, indC_min:indC_max) = boundary_local_new;
        end
    end
    myplot(boundary_new);

    %% After adding more gb and get ID_updated, match ID_updated to ID_iE, update ID and grain info
    % ID_updated: modified from the added grain boundaries
    % ID: the target ID map at this iE

    % umPerDp = 1;    % micron per data point in EBSD data
    % data_in.umPerDp = umPerDp;
    data_in.symmetry = 'hcp';
    data_in.ID_target = ID;
    data_in.ID_temp = ID_updated;
    data_in.x = x;
    data_in.y = y;
    data_in.phi1 = phi1;
    data_in.phi = phi;
    data_in.phi2 = phi2;
    data_in.gID = gID;
    data_in.gPhi1 = gPhi1;
    data_in.gPhi = gPhi;
    data_in.gPhi2 = gPhi2;

    data_out = match_ID_map_and_summarize(data_in);

    ID = data_out.ID;   % the modified ID map
    gID = data_out.gID;
    gPhi1 = data_out.gPhi1;
    gPhi = data_out.gPhi;
    gPhi2 = data_out.gPhi2;
    gCenterX = data_out.gCenterX;
    gCenterY = data_out.gCenterY;
    gNNeighbors = data_out.gNNeighbors;
    gDiameter = data_out.gDiameter;
    gArea = data_out.gArea;
    gEdge = data_out.gEdge;
    gNeighbors = data_out.gNeighbors;

    save_dir_2 = fullfile(save_dir, 'step-2');
    mkdir(save_dir_2);

    save(fullfile(save_dir_2, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']),'euler_aligned_to_sample','ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

    close all;

    disp(['finished twin grain file at iE=',num2str(iE)]);

end

%% copy to analysis folder
for iE = 0:iE_max
    try
        copyfile(fullfile(save_dir_2, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), ...
            fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), 'f');
    end
end

%% Part-3: Affine transform ID map. Link grains, and modify linked grains to the same ID#.
% Generate geotrans/tform information at iE = 1 to iE_max, without showing results
% Link ids (save in 'tbl') at differnt iEs after loading the saved geotrans/tform.
%   So, later, we can have grain files output again with different IDs.

% load updated variable 'grain_pair'
run(fullfile(p_setting,f_setting));

% reference, iE=0
d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));

x_0 = d.x;
y_0 = d.y;
ID_0 = d.ID;
boundary_0 = find_one_boundary_from_ID_matrix(ID_0);
gID_0 = d.gID;
gPhi1_0 = d.gPhi1;
gPhi_0 = d.gPhi;
gPhi2_0 = d.gPhi2;
gCenterX_0 = d.gCenterX;
gCenterY_0 = d.gCenterY;
gEdge_0 = d.gEdge;

for iE = 1:iE_max
    iB = iE + 1;
    disp(['iE=',num2str(iE)]);
    % load the modified EBSD data at iE
    load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'x','y','ID','gID','gEdge','gCenterX','gCenterY');

    try
        control_IDs = grain_pair([1,iB],:);
    catch
        control_IDs = [];
    end
    [tform, return_state] = align_ID_maps(ID_0, ID, control_IDs);
    tforms{iB} = tform;     % save the tform

    ID_0_to_iE = interp_data(x_0,y_0,ID_0, x,y,tform, 'interp', 'nearest');
    boundary_0_to_iE = find_boundary_from_ID_matrix(ID_0_to_iE);

    [ID_new, id_link_additional, id_link] = hungarian_assign_ID_map(ID_0_to_iE, ID);
    % remove 0s
    ind_r = find(id_link(:,1)==0 | id_link(:,2)==0);
    id_link(ind_r,:) = [];
    inds = isnan(ID_0_to_iE)|(ID_new==0); % old ID is nan, or matched ID=0
    ID_new(inds) = nan;

    % plot fine tformed map
    myplot(x,y, ID_new, boundary_0_to_iE);
    title(['iE=0 transformed with matched ID at iE=',num2str(iE)],'fontweight','normal');
    set(gca,'fontsize',18);
    print(fullfile(save_dir, ['B tform matched ID map iE=',num2str(iE),'.tif']),'-r300','-dtiff');
    close;

    tbl_iE = array2table(id_link);
    tbl_iE.Properties.VariableNames = {'iE_0',['iE_',num2str(iE)]};

    if iE==1
        tbl = tbl_iE;
    else
        tbl = innerjoin(tbl,tbl_iE);
    end
end

save(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms','tbl');

%% part-4, modify linked ID numbers to be the same throughout iEs
load(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms','tbl');

maxID = 0;
for iE = 0:iE_max
    f = matfile(fullfile(save_dir,[sample_name,'_parent_grain_file_iE_',num2str(iE)]));

    gID = f.gID;
    maxID = max(maxID, max(gID));
end
id_to_add = 10^(ceil(log10(maxID)));    % round to 1000 etc for adjustment.

save_dir_4 = fullfile(save_dir, 'step-4');
mkdir(save_dir_4);

for iE = 0:iE_max
    clear euler_aligned_to_sample;
    load(fullfile(save_dir,[sample_name,'_parent_grain_file_iE_',num2str(iE)]));

    if ~exist('euler_aligned_to_sample','var')
        warning('suggest to align_euler_to_sample first');
        euler_aligned_to_sample = 0;
    end

    % adjust, add a big enough value to prevent old unmodified grains having the same id as new modified grains
    ID(ID>0) = ID(ID>0) + id_to_add;
    gID(gID>0) = gID(gID>0) + id_to_add;
    gNeighbors(gNeighbors>0) = gNeighbors(gNeighbors>0) + id_to_add;    % add

    id_link = tbl.Variables;
    g_list_target = id_link(:,1);
    g_list_source = id_link(:,iE+1) + id_to_add;  % Temporarily change the gid in the source gID list.  Remember to use 'iE+1'
    for ii = 1:length(g_list_target)
        id_source = g_list_source(ii);
        id_target = g_list_target(ii);
        % modify
        ID(ID==id_source) = id_target;

        ind = find(gID_0==id_target);
        eulerd_ref = [gPhi1_0(ind), gPhi_0(ind), gPhi2_0(ind)];
        ind = find(gID==id_source);
        eulerd = [gPhi1(ind), gPhi(ind), gPhi2(ind)];

        % We can modify, or not modify. It won't affect the parient orientation determination in the twin identification.
        eulerd = find_closest_orientation_hcp(eulerd, eulerd_ref);      % ==> Find Closest Orientation

        gPhi1(gID==id_source) = eulerd(1);
        gPhi(gID==id_source) = eulerd(2);
        gPhi2(gID==id_source) = eulerd(3);

        gID(gID==id_source) = id_target;
        gNeighbors(gNeighbors==id_source) = id_target;
    end

    boundary = find_boundary_from_ID_matrix(ID);
    myplot(x,y,mod(ID,10),boundary);
    clim = caxis;
    %     caxis([0,maxID]);
    set(gca,'fontsize',18);
    title(['iE=',num2str(iE),' matched ID last digit'],'fontweight','normal');
    print(fullfile(save_dir, ['C matched ID map iE=',num2str(iE),'.tif']),'-r300','-dtiff');
    close;

    save(fullfile(save_dir_4, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'euler_aligned_to_sample', ...
        'ID','phi1','phi','phi2','x','y',...
        'gID','gPhi1','gPhi','gPhi2','gCenterX','gCenterY','gNNeighbors','gDiameter','gArea','gEdge','gNeighbors');

end

%%
for iE = 0:iE_max
    copyfile(fullfile(save_dir_4, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'f');
end


