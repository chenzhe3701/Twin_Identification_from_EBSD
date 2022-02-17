
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

%% (2.4) Merge, then draw mask to divide. copy after step-1

% run the .m file as setting file to load variables
run(fullfile(p_setting,f_setting));
assert(exist('ID_list','var')==1);
assert(exist('ID_merge_list','var')==1);

% target [iE,iD] list
iE_ID_list = [2,57;
    3,62];
iN_target = [];

iE = 0;
%% can run cell for each iE

iE
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

% So, here, the function is to :
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
    
    % only process the target [iE,ID]
    if ~ismember([iE,ID_current], iE_ID_list, 'rows') || (~isempty(iN_target) && iN~=iN_target)
        iN = iN + 1;
        continue;
    else
        mask_cell{iB}{iN} = [];
    end
    
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
    if (length(mask_cell{iB})>=iN) && ~isempty(mask_cell{iB}{iN})
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
            label_map_with_ID(x_local,y_local,ID_local, gcf, ID_current, 'r');
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
            caxism([5, 100]);
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


iE = iE + 1

