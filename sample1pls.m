%% Behaviour PLS Analysis with FEP data using civetsurf (https://github.com/katielavigne/civetsurf)
% Description:          Partial least squares assessing patterns of correlation between 
%                       cortical thickness (X) and cognition/demographics/symptoms (Y)
% Research Question:    What are the general and specific relationships between cognition and
%                       cortical thickness in FEP?
% Author: Katie Lavigne (katiemlavigne@gmail.com)

%% SETUP - Full Sample

% Data Preparation: prepare CIVET-processed structural MRI data & glimfile for analysis
data = dataprep(); % create data structure (no filter)

% Parcellation: create parcellation for vertex data
[data.parc] = parcellate(data.Y.smooth20mm, data.avsurf, data.mask, 'dkt', pwd); % parcellate

% Create domain-specific scores
data.cog = [data.glimfile.zVerbMem, data.glimfile.zVisMem, data.glimfile.zWorkMem, data.glimfile.zProcSpeed, data.glimfile.zExecFunc, data.glimfile.zAtt];
[data.cog_residmodel, ~, ~, data.cogresid] = ssregress(data.cog, data.glimfile, {'zGCI', 'Battery'});
newvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
for i = 1:size(newvars,2)
    data.glimfile.(newvars{i}) = data.cogresid(:,i);
end

data.gfields = fieldnames(data.glimfile);
save civetsurf_parc_fullsample.mat data % save output

%% VERTEX PLS
analysis_list = {'GenCT-GenCog', 'GenCT-SpecCog', 'SpecCT-GenCog', 'SpecCT-SpecCog'};

for i = 1:size(analysis_list,2)
    load civetsurf_parc_fullsample.mat
    % Convert Group to categorical
    data.group = categorical(data.glimfile.Group);
    data.gnames = categories(data.group);
    % Create group stacked data
    for g = 1:size(data.gnames,1)
        switch analysis_list{i}
            case 'GenCT-GenCog'
                covars = {'Battery'};
                behvars = {'zVerbMem', 'zVisMem', 'zWorkMem','zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education'};
%                 behvars = {'zVerbMem', 'zVisMem', 'zWorkMem','zProcSpeed', 'zExecFunc', 'zAtt'};
            case 'GenCT-SpecCog'
                covars = {'Battery'};
                behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education'};
%                 behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
            case 'SpecCT-GenCog'
                covars = {'Battery', 'meanCorticalMeasure20mm'};
                behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education'};
%                 behvars = {'zVerbMem', 'zVisMem', 'zWorkMem', 'zProcSpeed', 'zExecFunc', 'zAtt'};
            case 'SpecCT-SpecCog'
                covars = {'Battery', 'meanCorticalMeasure20mm'};
                behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec', 'Full_Scale_IQ', 'Age_at_Scan', 'Gender', 'Handed__ness', 'Years_of_Education'};
%                 behvars = {'zVerbMem_spec', 'zVisMem_spec', 'zWorkMem_spec', 'zProcSpeed_spec', 'zExecFunc_spec', 'zAtt_spec'};
        end
        behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'IQ', 'Age', 'Sex', 'Handedness', 'Education'};
%         behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention'};
        % Regress out covariates
        rows = data.group==data.gnames{g};
        [data.(['residmodel_' data.gnames{g}]), ~, ~, data.(['resid_' data.gnames{g}])] = ssregress(data.Y.smooth20mm(rows,:), data.glimfile(rows,:), covars);
        data.stacked_data{g,1} = data.(['resid_' data.gnames{g}]);
    end
    
    % Run PLS
    PLS = pls(data.stacked_data, data, behvars, behdesc); % run pls
    % Save PLS analysis
    save(fullfile('pls', ['civetsurf_pls_' analysis_list{i} '.mat'])) % save results output
    movefile('pls', ['pls_vertex_' analysis_list{i}]) % rename directory
    close all % close figures
end