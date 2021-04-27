%% Behaviour PLS Analysis with FEP data using civetsurf (https://github.com/katielavigne/civetsurf)
% Description:          Partial least squares assessing patterns of correlation between 
%                       cortical thickness (X) and cognition/demographics/symptoms (Y)
% Research Question:    What are the general and specific relationships between cognition and
%                       cortical thickness in FEP?
% Analysis Workflow
%   1. Setup
%   2. Behaviour PLS analyses covarying for battery
%       a. Global CT - General Cognition
%       b. Global CT - Specific Cognition (covary general cognition)
%       c. Regional CT - General Cognition (covary mean thickness)
%       d. Regional CT - Specific Cognition (covary both)

% Author: Katie Lavigne (katiemlavigne@gmail.com)

%% SETUP

% Data Preparation: prepare CIVET-processed structural MRI data & glimfile for analysis
data = dataprep(); % create data structure (filter: Group.FEP)

% Parcellation: create parcellation for vertex data.
[data.parc] = parcellate(data.Y.smooth20mm, data.avsurf, data.mask, 'dkt', pwd); % parcellate

% Create domain-specific scores
data.cog = [data.glimfile.zVerbMem, data.glimfile.zVisMem, data.glimfile.zWorkMem, data.glimfile.zProcSpeed, data.glimfile.zExecFunc, data.glimfile.zAtt];
[data.cog_residmodel, ~, ~, data.cogresid] = ssregress(data.cog, data.glimfile, {'zGCI', 'Battery'});
newvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
for i = 1:size(newvars,2)
    data.glimfile.(newvars{i}) = data.cogresid(:,i);
end

save civetsurf_parc.mat data % save output

%% PLS - Sample 1 - ROI

analysis_list = {'GenCT-GenCog', 'GenCT-SpecCog', 'SpecCT-GenCog', 'SpecCT-SpecCog'};

for i = 1:size(analysis_list,2)
    load civetsurf_parc.mat
    switch analysis_list{i}
        case 'GenCT-GenCog'
            covars = {'Battery'};
            behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
        case 'GenCT-SpecCog'
            covars = {'Battery'};
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
        case 'SpecCT-GenCog'
            covars = {'Battery', 'meanCorticalMeasure20mm'};
            behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
        case 'SpecCT-SpecCog'
            covars = {'Battery', 'meanCorticalMeasure20mm'};
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
    end
    behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'IQ', 'Age', 'Sex', 'Handedness', 'Education', 'SAPS Total', 'SANS Total', 'CPZ Adherence'};
    % Regress out covariates
    [data.residmodel, ~, ~, data.resid] = ssregress(data.parc.ROIs, data.glimfile, covars); 
    % Run PLS
    PLS = pls(data.resid, data.parc.pinfo, data.mask, data.avsurf, data.glimfile, behvars, behdesc); % run pls
    % Save PLS analysis
    save(fullfile('pls', ['civetsurf_pls_' analysis_list{i} '.mat'])) % save results output
    movefile('pls', ['pls_' analysis_list{i}]) % rename directory
    close all % close figures
end

%% VERTEX PLS
analysis_list = {'GenCT-GenCog', 'GenCT-SpecCog', 'SpecCT-GenCog', 'SpecCT-SpecCog'};

for i = 1:size(analysis_list,2)
    load civetsurf_parc.mat
    switch analysis_list{i}
        case 'GenCT-GenCog'
            covars = {'Battery'};
            behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
        case 'GenCT-SpecCog'
            covars = {'Battery'};
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
        case 'SpecCT-GenCog'
            covars = {'Battery', 'meanCorticalMeasure20mm'};
            behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
        case 'SpecCT-SpecCog'
            covars = {'Battery', 'meanCorticalMeasure20mm'};
            behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age', 'Gender', 'Handedness', 'Years_of_Education', 'SAPS_Total', 'SANS_Total_No_Att', 'CPZ_with_Adherence'};
    end
    behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'IQ', 'Age', 'Sex', 'Handedness', 'Education', 'SAPS Total', 'SANS Total', 'CPZ_Adherence'};
    % Regress out covariates
    [data.residmodel, ~, ~, data.resid] = ssregress(data.Y.smooth20mm, data.glimfile, covars); 
    % Run PLS
    PLS = pls(data.resid, data.parc.pinfo, data.mask, data.avsurf, data.glimfile, behvars, behdesc); % run pls
    % Save PLS analysis
    save(fullfile('pls', ['civetsurf_pls_' analysis_list{i} '.mat'])) % save results output
    movefile('pls', ['pls_vertex_' analysis_list{i}]) % rename directory
    close all % close figures
end