
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

%% part-5: find out variants
run(fullfile(p_setting,f_setting));
assert(exist('min_gs','var')==1);

po_tolerance_angle = 10; % if child grain has misorientation < po_tolerance_angle with undeformed parent grain, it is considered as having parent orientation
twin_tolerance_angle = 10;  % if child grain has misorientation < twin_tolerance_angle to a potential twin variant, it is identified as that twin variant
for iE = 0:iE_max
    iB = iE + 1;

    d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_0.mat']));
    gID_0 = d.gID;
    gPhi1_0 = d.gPhi1;
    gPhi_0 = d.gPhi;
    gPhi2_0 = d.gPhi2;
    ID_0 = d.ID;
    boundary_0 = find_boundary_from_ID_matrix(ID_0);

    d = load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']));
    x = d.x;
    y = d.y;
    gID_p = d.gID;
    gPhi1_p = d.gPhi1;
    gPhi_p = d.gPhi;
    gPhi2_p = d.gPhi2;
    ID_p = d.ID;
    boundary_p = find_boundary_from_ID_matrix(ID_p);

    % data where twins (children) are individially labeled with IDs
    d = load(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']));
    gID_c = d.gID;
    gPhi1_c = d.gPhi1;
    gPhi_c = d.gPhi;
    gPhi2_c = d.gPhi2;
    ID_c = d.ID;
    phi1_c = d.phi1;
    phi_c = d.phi;
    phi2_c = d.phi2;

    ID_variant_grain_wise = zeros(size(ID_p));
    ID_variant_point_wise = zeros(size(ID_p));
    Misorientation_point_wise = zeros(size(ID_p));

    for ii = 1: length(gID_p)
        id_p = gID_p(ii);
        disp(['current ID = ',num2str(id_p)]);

        % Find the parent orientation at iE = 0
        ind_0 = find(gID_0 == id_p);
        euler_0 = [gPhi1_0(ind_0), gPhi_0(ind_0), gPhi2_0(ind_0)];

        if ~isempty(ind_0)
            % find out all the id numbers on ID_c overlap with the parent grain on ID_p
            id_c = unique(ID_c(ID_p == id_p));

            % Additionally, the child should have > 70% of its pixels
            % overlap with the potential parent grain
            % ==> 2021-11-19 I don't think we need this anymore

            % Modify all the child orientation to crystallographically equivalent one that is closest to euler_0
            misorientation = [];
            col_ID = [];
            col_euler = [];
            for jj = 1:length(id_c)
                id = id_c(jj);
                ind = (gID_c == id);
                euler_id = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)];

                % =================> Important, here modified the grain info of the child grain
                euler_id = find_closest_orientation_hcp(euler_id, euler_0);
                gPhi1_c(ind) = euler_id(1);
                gPhi_c(ind) = euler_id(2);
                gPhi2_c(ind) = euler_id(3);

                misorientation(jj,1) = calculate_misorientation_euler_d(euler_0, euler_id, 'hcp');
                col_ID = [col_ID; id];
                col_euler = [col_euler; euler_id];
            end
            tbl = table(col_ID, col_euler, misorientation); % for debugging

            % (1) Find out which 'children' is actually the 'parent', by checking the misorientation
            % Note: multiple 'grains' can have the parent orientation
            inds = find(misorientation < po_tolerance_angle);
            if ~isempty(inds)
                % Use the average of all the 'child grains' with parent orientation.
                id_po = id_c(inds);
                inds_po = find(ismember(gID_c,id_po));
                euler_po = calculate_average_dominant_euler_hcp([gPhi1_c(inds_po), gPhi_c(inds_po), gPhi2_c(inds_po)]);

                % remove the 'child grains' with the 'parent orientation', leave only the true 'twin children'
                id_c(inds) = [];
            else
                % the child maybe fully twinned. Just use euler_0 as euler_po
                euler_po = euler_0;
                warning('The grain might have completely twinned');
                str = sprintf('iE = %d, ii = %d, ID_parent = %d \n', iE, ii, id_p);
                disp(str);
            end

            % (2) Find out if the remaining child grain is a twin
            for jj = 1:length(id_c)
                id = id_c(jj);  % twin grains id
                ind = (gID_c == id);
                euler_id = [gPhi1_c(ind), gPhi_c(ind), gPhi2_c(ind)];

                misorientation = [];
                %  if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                %       hcp_cell('euler',euler_id,'material','Mg','plotPlane',0,'plotBurgers',0,'plotTrace',0);
                %  end

                % Determine if this child grain is a twin variant:
                % Compare [euler of this child] vs [euler of the iTwin system of the parent orientation]
                for kk = 1:6
                    euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                    misorientation(kk) = calculate_misorientation_euler_d(euler_id, euler_kk, 'HCP');
                    % if (demoTF==1)&&(jj==1)&&(ii==ii_demo)
                    % 	hcp_cell('euler',euler_po,'material','Mg','ss',24+kk);
                    % 	hcp_cell('euler',euler_kk,'material','Mg','ss',24+kk,'plotPlane',0,'plotBurgers',0,'plotTrace',0);
                    % end
                end
                [min_val, iVariant_child] = min(abs(misorientation));    % find out

                % ==============> The child grain may be a twin area containing multiple variants. Assume the child orientation represents at least one true twin orientation
                % If small enough, the child grain should be a twin. Do point-wise analysis
                if min_val < twin_tolerance_angle && sum(ind) >= min_gs
                    ID_variant_grain_wise(ID_c == id) = iVariant_child; % grain-wise variant map

                    ind_list = find(ID_c==id);
                    for kk = 1:length(ind_list)
                        ind = ind_list(kk);
                        euler_c = [phi1_c(ind), phi_c(ind), phi2_c(ind)];

                        misorientation = [];
                        % compare [euler of this pixel]  vs  [euler of the iTwin system of the parent orientation]
                        for kk = 1:6
                            euler_kk = euler_by_twin(euler_po, kk, 'Mg');    % calculate the twin euler angle for twin system kk
                            misorientation(kk) = calculate_misorientation_euler_d(euler_c, euler_kk, 'HCP');
                        end
                        [miso, iVariant] = min(abs(misorientation));
                        Misorientation_point_wise(ind) = miso;
                        ID_variant_point_wise(ind) = iVariant;
                    end
                elseif sum(ind) < min_gs
                    warning('child grain smaller than min_gs, not considered as a twin');
                    str = sprintf('iE = %d, ii = %d, ID_parent = %d, jj = %d, ID_twin = %d \n', iE, ii, id_p, jj, id);
                    disp(str);
                else
                    warning('Twin grain misorientation with parent orientation > 10 deg, rejected as a variant:');
                    str = sprintf('iE = %d, ii = %d, ID_parent = %d, jj = %d, ID_twin = %d \n', iE, ii, id_p, jj, id);
                    disp(str);
                end
            end
        end
    end
    variant_grain_wise{iB} = ID_variant_grain_wise;
    variant_point_wise{iB} = ID_variant_point_wise;

    % After processing this iE, Update the modified CHILD grain data and save
    save_dir_5 = fullfile(save_dir, 'step-5');
    mkdir(save_dir_5);

    copyfile(fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir_5, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), 'f');
    gPhi1 = gPhi1_c;
    gPhi = gPhi_c;
    gPhi2 = gPhi2_c;
    save(fullfile(save_dir_5, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), 'gPhi1','gPhi','gPhi2', '-append');


    myplot(x,y, ID_variant_grain_wise, boundary_p); caxis([0 6]);
    make_variant_map_background('tickOff',false);
    set(gca,'fontsize',18);
    title(['iE=',num2str(iE)],'fontweight','normal');
    print(fullfile(save_dir,['variant_grain_wise_iE=',num2str(iE),'.tif']),'-dtiff');
    close;

    myplot(x,y, ID_variant_point_wise, boundary_p); caxis([0 6]);
    make_variant_map_background('tickOff',false);
    set(gca,'fontsize',18);
    title(['iE=',num2str(iE)],'fontweight','normal');
    print(fullfile(save_dir,['variant_pt_wise_iE=',num2str(iE),'.tif']),'-dtiff');
    close;

end

save(fullfile(save_dir,'variant_maps.mat'),'variant_grain_wise','variant_point_wise');

% copy the updated file back to the main folder
for iE = 0:iE_max
    copyfile(fullfile(save_dir_5, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), ...
        fullfile(save_dir, [sample_name,'_grain_file_iE_',num2str(iE),'.mat']), 'f');
end


%% part-6: summarize twin pct
load(fullfile(save_dir,'variant_maps.mat'),'variant_grain_wise','variant_point_wise');
nr = 3;
nc = 3;
twinPct = zeros(nr*nc,iE_max+1);
for iE = 0:iE_max
    iB = iE + 1;

    load(fullfile(save_dir, [sample_name,'_parent_grain_file_iE_',num2str(iE),'.mat']), 'ID');

    variant_map = variant_point_wise{iB};
    if iB==1
        [nR,nC] = size(variant_map);
    end
    for ir=1:nr
        for ic = 1:nc
            sub_variant_map = variant_map([1:floor(nR/nr)] + floor(nR/nr)*(ir-1), [1:floor(nC/nc)] + floor(nC/nc)*(ic-1));
            sub_ID_map = ID([1:floor(nR/nr)] + floor(nR/nr)*(ir-1), [1:floor(nC/nc)] + floor(nC/nc)*(ic-1)); % just for counting number of pixels

            ii = iE + 1;
            twinPct((ir-1)*nc+ic,ii) = sum(sub_variant_map(:)>0)/sum(sub_ID_map(:)>0);
        end
    end
end

tAvg = mean(twinPct);
tStd = std(twinPct);

save(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd');

%% part-7, calculate EBSD estimated strain
load(fullfile(save_dir,'geotrans_and_id_link.mat'),'tforms');
for iE = 0:iE_max
    iB = iE+1;
    if ~isempty(tforms{iB})
        epsilon = tform_to_epsilon(tforms{iB});
        strain_ebsd(iB) = round(epsilon(1), 4);
    else
        strain_ebsd(iB) = 0;
    end
end
str = [sprintf('strain_ebsd = ['), sprintf('%.4f, ',strain_ebsd(1:4)), newline, ...
    sprintf('%.4f, ',strain_ebsd(5:8)), newline, ...
    sprintf('%.4f, ',strain_ebsd(9:11)), newline, ...
    sprintf('%.4f, ',strain_ebsd(12:13)), sprintf('%.4f];',strain_ebsd(14))];
disp(str);

% run the .m file as setting file to load variables
run(fullfile(p_setting,f_setting));
assert(exist('strain_sg','var')==1);

save(fullfile(save_dir, 'twin_pct.mat'), 'strain_ebsd', 'strain_sg', '-append');
%% plot
load(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd', 'strain_ebsd', 'strain_sg');

% run the .m file as setting file to load variables
run(fullfile(p_setting,f_setting));
assert(exist('inds_half_cycle','var')==1);

N = length(inds_half_cycle);
colors = parula(N+1);
% [1] using strain gage strain
close all;
figure; hold on;
clear inds;
for ii = 1:length(inds_half_cycle)
    if ii == 1
        inds{ii} = 1:inds_half_cycle(ii);
    else
        inds{ii} = inds_half_cycle(ii-1):inds_half_cycle(ii);
    end
    errorbar(strain_sg(inds{ii}), 100*tAvg(inds{ii}), 100*tStd(inds{ii}), '.-', 'color',colors(ii,:), 'linewidth',1.5,'markersize',24);
end
ymin = round(min((tAvg-tStd)*100),0) - 2;
ymax = round((max(tAvg+tStd)*100),-1) + 10;
set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[-0.04, 0.005],'ylim',[ymin ymax],'fontsize',18,'fontweight','normal');
xlabel('Strain from strain gage');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_sg.tiff'),'-dtiff');

% [2] using (fine transformed) ebsd estimated strain
figure; hold on;
for ii = 1:length(inds_half_cycle)
    errorbar(strain_ebsd(inds{ii}), 100*tAvg(inds{ii}), 100*tStd(inds{ii}), '.-', 'color',colors(ii,:), 'linewidth',1.5,'markersize',24);
end
xmin = round(min(strain_ebsd),2)-0.01;
xmax = round(max(strain_ebsd),2)+0.01;
set(gca,'xdir','normal','linewidth',1.5);
set(gca,'xlim',[xmin, xmax],'ylim',[ymin, ymax],'fontsize',18,'fontweight','normal');
xlabel('Strain from ebsd estimate');
ylabel('Twin Area Percent (%)');
print(fullfile(save_dir,'twin_pct_vs_ebsd_strain.tiff'),'-dtiff');


tbl = array2table([(0:iE_max)', strain_sg(:), strain_ebsd(:), 100*tAvg(:), 100*tStd(:)]);
tbl.Properties.VariableNames = {'iE','strain_sg','strain_ebsd','twinPct %','twinStd %'};
disp(tbl);
figure;
uitable('Data',tbl{:,:},'ColumnName',tbl.Properties.VariableNames,...
    'RowName',tbl.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
print(fullfile(save_dir,'twin pct table.tiff'),'-dtiff');

save(fullfile(save_dir, 'twin_pct.mat'), 'twinPct', 'tAvg', 'tStd', 'tbl');

close all;
