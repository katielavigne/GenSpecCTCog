%% Behaviour PLS Analysis with FEP data using civetsurf (https://github.com/katielavigne/civetsurf)
% Description:          Partial least squares assessing patterns of correlation between 
%                       cortical thickness (X) and cognition/demographics/symptoms (Y)
% Research Question:    What are the general and specific relationships between cognition and
%                       cortical thickness in FEP?
% Author: Katie Lavigne (katiemlavigne@gmail.com)

%% SETUP - SAMPLE 1

% Data Preparation: prepare CIVET-processed structural MRI data & glimfile for analysis
data1 = dataprep(); % create data structure (filter: Group.Patient)

% Create domain-specific scores
data1.cog = [data1.glimfile.zVerbMem, data1.glimfile.zVisMem, data1.glimfile.zWorkMem, data1.glimfile.zProcSpeed, data1.glimfile.zExecFunc, data1.glimfile.zAtt];
[data1.cog_residmodel, ~, ~, data1.cogresid] = ssregress(data1.cog, data1.glimfile, {'zGCI'});
newvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
for i = 1:size(newvars,2)
    data1.glimfile.(newvars{i}) = data1.cogresid(:,i);
end

data1.gfields = fieldnames(data1.glimfile);

%% SETUP - Sample 2

% Data Preparation: prepare CIVET-processed structural MRI data & glimfile for analysis.
data2 = dataprep(); % create data structure (filter: group.Patient)

% Create domain-specific scores
data2.cog = [data2.glimfile.cog_verb_mem_z, data2.glimfile.cog_vis_mem_z, data2.glimfile.cog_work_mem_z, data2.glimfile.cog_speed_proc_z, data2.glimfile.cog_exec_func_z, data2.glimfile.cog_vis_att_z];
data2.glimfile.zGCI = mean(data2.cog,2);
[data2.cog_residmodel, ~, ~, data2.cogresid] = ssregress(data2.cog, data2.glimfile, {'zGCI'});
newvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
for i = 1:size(newvars,2)
    data2.glimfile.(newvars{i}) = data2.cogresid(:,i);
end

data2.gfields = fieldnames(data2.glimfile);

%% COMBAT HARMONIZATION

Y1 = data1.Y.smooth20mm';
Y2 = data2.Y.smooth20mm';
Y = [Y1 Y2];

b1 = zeros(1,size(Y1,2));
b2 = data2.glimfile.Pre_Post_MRI__upgrade_date__20180917_';
batch = [b1 b2] + 1;

a1 = data1.glimfile.Age_at_Scan;
a2 = data2.glimfile.dem_age;
age = [a1; a2];

gender = [data1.glimfile.Gender; data2.glimfile.sex];
gender = dummyvar(gender);

group = [data1.glimfile.Group; data2.glimfile.group];
group = categorical(group);
group = dummyvar(group);

mod = [age gender(:,2)];

data_harmonized = combat(Y, batch, mod, 1);
data1.harmonized = data_harmonized(:,1:97)';
data2.harmonized = data_harmonized(:,98:end)';

data1.glimfile.mean_thickness20mm_2_1_1_harmonized = mean(data1.harmonized,2);
data2.glimfile.mean_thickness20mm_2_1_1_harmonized = mean(data2.harmonized,2);

% Parcellation: create parcellation for vertex data
[data1.parc] = parcellate(data1.harmonized, data1.avsurf, data1.mask, 'dkt', pwd); % parcellate
[data2.parc] = parcellate(data2.harmonized, data2.avsurf, data2.mask, 'dkt', pwd); % parcellate

save civetsurf_sample1.mat data1 % save output
save civetsurf_sample2.mat data2 % save output

% FIX SCAN NOTES AND OTHER STRINGS WITH RETURNS
for i = 1:size(data2.glimfile,2)
    if iscell(data2.glimfile.(data2.gfields{i}))
        data2.glimfile.(data2.gfields{i}) = regexprep(data2.glimfile.(data2.gfields{i}), '\n', ';');
    end
end


%% SAMPLE 1 VERTEX PLS
analysis_list = {'GenCT-GenCog', 'GenCT-SpecCog', 'SpecCT-GenCog', 'SpecCT-SpecCog'};

for i = 1:size(analysis_list,2)
    load civetsurf_sample1.mat
    switch analysis_list{i}
        case 'GenCT-GenCog'
            covars = {'Battery'};
            behvars = {'zVerbMem', 'zVisMem', 'zWorkMem','zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education', 'SANS_Total_No_Att', 'SAPS_Total', 'CPZ_with_Adherence_at_Scan'};
%                 behvars = {'zVerbMem', 'zVisMem', 'zWorkMem','zProcSpeed', 'zExecFunc', 'zAtt'};
        case 'GenCT-SpecCog'
            covars = {'Battery'};
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education', 'SANS_Total_No_Att', 'SAPS_Total', 'CPZ_with_Adherence_at_Scan'};
%                 behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
        case 'SpecCT-GenCog'
            covars = {'Battery', 'mean_thickness20mm_2_1_1_harmonized'};
            behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education', 'SANS_Total_No_Att', 'SAPS_Total', 'CPZ_with_Adherence_at_Scan'};
%                 behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt'};
        case 'SpecCT-SpecCog'
            covars = {'Battery', 'mean_thickness20mm_2_1_1_harmonized'};
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education', 'SANS_Total_No_Att', 'SAPS_Total', 'CPZ_with_Adherence_at_Scan'};
%                 behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
    end
    behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'IQ', 'Age', 'Sex', 'Handedness', 'Education', 'SANS Total', 'SAPS Total', 'CPZ Adherence'};
%         behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention'};
    % Regress out covariates
    [data1.residmodel, ~, ~, data1.resid] = ssregress(data1.harmonized, data1.glimfile, covars);    
    % Run PLS
    PLS = pls(data1.resid, data1, behvars, behdesc); % run pls
    % Save PLS analysis
    save(fullfile('pls', ['civetsurf_pls_' analysis_list{i} '.mat'])) % save results output
    movefile('pls', ['pls1_vertex_' analysis_list{i}]) % rename directory
    close all % close figures
end

%% SAMPLE 2 VERTEX PLS
analysis_list = {'GenCT-GenCog', 'GenCT-SpecCog', 'SpecCT-GenCog', 'SpecCT-SpecCog'};

for i = 1:size(analysis_list,2)
    load civetsurf_sample2.mat
    switch analysis_list{i}
        case 'GenCT-GenCog'
            covars = '';
            behvars = {'cog_verb_mem_z', 'cog_vis_mem_z', 'cog_work_mem_z', 'cog_speed_proc_z', 'cog_exec_func_z', 'cog_vis_att_z', 'wasi_full_iq_TP1', 'dem_age', 'sex', 'hand_bin_TP1', 'dem_years_educ', 'sans_total_no_attention', 'saps_total', 'CPZ_Adherence'};
%                 behvars = {'cog_verb_mem_z', 'cog_vis_mem_z', 'cog_work_mem_z', 'cog_speed_proc_z', 'cog_exec_func_z', 'cog_vis_att_z'};
        case 'GenCT-SpecCog'
            covars = '';
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'wasi_full_iq_TP1', 'dem_age', 'sex', 'hand_bin_TP1', 'dem_years_educ', 'sans_total_no_attention', 'saps_total', 'CPZ_Adherence'};
%                 behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
        case 'SpecCT-GenCog'
            covars = {'mean_thickness20mm_2_1_1_harmonized'};
            behvars = {'cog_verb_mem_z', 'cog_vis_mem_z', 'cog_work_mem_z', 'cog_speed_proc_z', 'cog_exec_func_z', 'cog_vis_att_z', 'wasi_full_iq_TP1', 'dem_age', 'sex', 'hand_bin_TP1', 'dem_years_educ', 'sans_total_no_attention', 'saps_total', 'CPZ_Adherence'};
%                 behvars = {'cog_verb_mem_z', 'cog_vis_mem_z', 'cog_work_mem_z', 'cog_speed_proc_z', 'cog_exec_func_z', 'cog_vis_att_z'};
        case 'SpecCT-SpecCog'
            covars = {'mean_thickness20mm_2_1_1_harmonized'};
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'wasi_full_iq_TP1', 'dem_age', 'sex', 'hand_bin_TP1', 'dem_years_educ', 'sans_total_no_attention', 'saps_total', 'CPZ_Adherence'};
%                 behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
    end
    behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'IQ', 'Age', 'Sex', 'Handedness', 'Education', 'SANS Total', 'SAPS Total', 'CPZ Adherence'};
%         behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention'};
        % Regress out covariates
        if ~isempty(covars)
            [data2.residmodel, ~, ~, data2.resid] = ssregress(data2.harmonized, data2.glimfile, covars);
            PLS = pls(data2.resid, data2, behvars, behdesc); % run pls
        else
            PLS = pls(data2.harmonized, data2, behvars, behdesc); % run pls
        end
    % Save PLS analysis
    save(fullfile('pls', ['civetsurf_pls_' analysis_list{i} '.mat'])) % save results output
    movefile('pls', ['pls2_vertex_' analysis_list{i}]) % rename directory
    close all % close figures
end