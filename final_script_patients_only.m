%% Behaviour PLS Analysis with FEP data using civetsurf (https://github.com/katielavigne/civetsurf)
% Description:          Partial least squares assessing patterns of correlation between 
%                       cortical thickness (X) and cognition/demographics/symptoms (Y)
% Research Question:    What are the general and specific relationships between cortical thickness and
%                       cognition in FEP?
% Author: Katie Lavigne

%% SETUP - SAMPLE 1

% Data Preparation: prepare CIVET-processed structural MRI data & glimfile for analysis
data1 = dataprep(); % create data structure

% Create domain-specific scores
data1.cog = [data1.glimfile.z_vm, data1.glimfile.z_vism, data1.glimfile.z_wm, data1.glimfile.z_ps, data1.glimfile.z_ef, data1.glimfile.z_att];
[data1.cog_residmodel, ~, ~, data1.cogresid] = ssregress(data1.cog, data1.glimfile, {'z_gci'});
newvars = {'z_vm_spec', 'z_vism_spec', 'z_wm_spec', 'z_ps_spec', 'z_ef_spec', 'z_att_spec'};
for i = 1:size(newvars,2)
    data1.glimfile.(newvars{i}) = data1.cogresid(:,i);
end

data1.gfields = fieldnames(data1.glimfile);
movefile('Global cortical thickness measures.csv', 'Data 1 Global cortical thickness measures.csv')

%% SETUP - Sample 2

% Data Preparation: prepare CIVET-processed structural MRI data & glimfile for analysis.
data2 = dataprep(); % create data structure

% Create domain-specific scores
data2.cog = [data2.glimfile.z_vm, data2.glimfile.z_vism, data2.glimfile.z_wm, data2.glimfile.z_ps, data2.glimfile.z_ef, data2.glimfile.z_att];
data2.glimfile.z_gci= mean(data2.cog,2);
[data2.cog_residmodel, ~, ~, data2.cogresid] = ssregress(data2.cog, data2.glimfile, {'z_gci'});
newvars = {'z_vm_spec', 'z_vism_spec', 'z_wm_spec', 'z_ps_spec', 'z_ef_spec', 'z_att_spec'};
for i = 1:size(newvars,2)
    data2.glimfile.(newvars{i}) = data2.cogresid(:,i);
end

data2.gfields = fieldnames(data2.glimfile);
movefile('Global cortical thickness measures.csv', 'Data 2 Global cortical thickness measures.csv')

%% COMBAT HARMONIZATION

Y1 = data1.Y.smooth20mm';
Y2 = data2.Y.smooth20mm';
Y = [Y1 Y2];

b1 = zeros(1,size(Y1,2));
b2 = data2.glimfile.prepost';
batch = [b1 b2] + 1;

age = [data1.glimfile.age; data2.glimfile.age];

gender = [data1.glimfile.sex; data2.glimfile.sex];
gender = dummyvar(gender);

group = [data1.glimfile.group; data2.glimfile.group];
group = categorical(group);
group = dummyvar(group);

mod = [age gender(:,2)];

data_harmonized = combat(Y, batch, mod, 1);
data1.harmonized = data_harmonized(:,1:size(Y1,2))';
data2.harmonized = data_harmonized(:,size(Y1,2)+1:end)';

data1.glimfile.mean_thickness20mm_2_1_1_harmonized = mean(data1.harmonized,2);
data2.glimfile.mean_thickness20mm_2_1_1_harmonized = mean(data2.harmonized,2);

% Parcellation: create parcellation for vertex data
[data1.parc] = parcellate(data1.harmonized, data1.avsurf, data1.mask, 'dkt', pwd); % parcellate
[data2.parc] = parcellate(data2.harmonized, data2.avsurf, data2.mask, 'dkt', pwd); % parcellate

save civetsurf_sample1.mat data1 % save output
save civetsurf_sample2.mat data2 % save output
writetable(data1.glimfile, 'data1_glimfile.csv')
writetable(data2.glimfile, 'data2_glimfile.csv')

% Patients only
d1p = data1;
[d1p.filter, d1p.glimfile] = filter_data(d1p.glimfile, d1p.gfields);
d1p.harmonized = d1p.harmonized(d1p.prefilter.glimfile.group==2,:);
save civetsurf_sample1_patients.mat d1p

d2p = data2;
[d2p.filter, d2p.glimfile] = filter_data(d2p.glimfile, d2p.gfields);
d2p.harmonized = d2p.harmonized(d2p.prefilter.glimfile.group==2,:);
save civetsurf_sample2_patients.mat d2p
%% SAMPLE 1 VERTEX PLS

analysis_list = {'GenCT-Cog', 'SpecCT-Cog'}; 

for i = 1:size(analysis_list,2)
    load civetsurf_sample1_patients.mat
    switch analysis_list{i}
        case 'GenCT-Cog'
            covars = {'battery'};
            behvars = {'z_vm', 'z_vism', 'z_wm','z_ps', 'z_ef', 'z_att', 'iq', 'age', 'sex', 'handedness', 'education', 'sans', 'saps'};
        case 'SpecCT-Cog'
            covars = {'battery', 'mean_thickness20mm_2_1_1_harmonized'};
            behvars = {'z_vm', 'z_vism', 'z_wm','z_ps', 'z_ef', 'z_att', 'iq', 'age', 'sex', 'handedness', 'education', 'sans', 'saps'};
    end
    behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'IQ', 'Age', 'Sex', 'Handedness', 'Education', 'Negative Symptoms', 'Positive Symptoms'};
    % Regress out covariates
    [d1p.residmodel, ~, ~, d1p.resid] = ssregress(d1p.harmonized, d1p.glimfile, covars);    
    % Run PLS
    PLS = pls(d1p.resid, d1p, behvars, behdesc); % run pls
    % Save PLS analysis
    save(fullfile('pls', ['civetsurf_pls_' analysis_list{i} '.mat'])) % save results output
    movefile('pls', ['pls1_vertex_' analysis_list{i}]) % rename directory
    close all % close figures
end

%% SAMPLE 2 VERTEX PLS
analysis_list = {'GenCT-Cog', 'SpecCT-GenCog'};

for i = 1:size(analysis_list,2)
    load civetsurf_sample2_patients.mat
    switch analysis_list{i}
        case 'GenCT-Cog'
            covars = '';
            behvars = {'z_vm', 'z_vism', 'z_wm', 'z_ps', 'z_ef', 'z_att', 'iq', 'age', 'sex', 'handedness', 'education', 'sans', 'saps'};
        case 'SpecCT-Cog'
            covars = {'mean_thickness20mm_2_1_1_harmonized'};
            behvars = {'z_vm', 'z_vism', 'z_wm', 'z_ps', 'z_ef', 'z_att', 'iq', 'age', 'sex', 'handedness', 'education', 'sans', 'saps'};
    end
    behdesc = {'Verbal Memory', 'Visual Memory', 'Working Memory', 'Processing Speed', 'Executive Function', 'Attention', 'IQ', 'Age', 'Sex', 'Handedness', 'Education', 'Negative Symptoms', 'Positive Symptoms'};
        % Regress out covariates
        if ~isempty(covars)
            [d2p.residmodel, ~, ~, d2p.resid] = ssregress(d2p.harmonized, d2p.glimfile, covars);
            PLS = pls(d2p.resid, d2p, behvars, behdesc); % run pls
        else
            PLS = pls(d2p.harmonized, d2p, behvars, behdesc); % run pls
        end
    % Save PLS analysis
    save(fullfile('pls', ['civetsurf_pls_' analysis_list{i} '.mat'])) % save results output
    movefile('pls', ['pls2_vertex_' analysis_list{i}]) % rename directory
    close all % close figures
end

%% EXTRACT LV1 SpecCT-Cog sample 1 (ps = 0.052)

load pls1_vertex_SpecCT-Cog\civetsurf_pls_SpecCT-Cog.mat
mkdir pls/flipped

plsbar(PLS.result, 1, 1, behdesc, {}) % behavioural data
plssurf(d2p.parc.pinfo, d2p.mask, d2p.avsurf, PLS.result, 1, 1, [1.96, 2.58], 'Bootstrap Ratios') % brain data

% flipped
plsbar(PLS.result, 1, -1, behdesc, {}) % behavioural data
plssurf(d2p.parc.pinfo, d2p.mask, d2p.avsurf, PLS.result, 1, -1, [1.96, 2.58], 'Bootstrap Ratios') % brain data

movefile('pls/flipped/*', 'pls1_vertex_SpecCT-Cog\flipped\')
movefile('pls/LV*', 'pls1_vertex_SpecCT-Cog\')
rmdir('pls', 's')
close all

%% FIGURE 2 - GenCT-Cog

% PREP
load pls1_vertex_GenCT-Cog\civetsurf_pls_GenCT-Cog.mat
pls1 = PLS;
avsurf = data1.avsurf;
mask = data1.mask;
load pls2_vertex_GenCT-Cog\civetsurf_pls_GenCT-Cog.mat
pls2 = PLS;
clearvars -except pls1 pls2 avsurf mask behdesc

% SURFSTAT OVERLAY 
U1 = pls1.result.boot_result.compare_u(:,1); % comp1
U2 = pls2.result.boot_result.compare_u(:,1); % comp1
surfstat_overlay(U1,U2,mask,avsurf)
SurfStatColLim([1,4])
cmap = [0.8 0.8 0.8; cool(3)];
SurfStatColormap(cmap)
saveas(gcf, 'GenCT-Cog_overlay.fig')
saveas(gcf, 'GenCT-Cog_overlay.png')
close all

% BEH OVERLAY
pls_overlay(pls1, 1, pls2, 1, behdesc, [1, 1])
saveas(gcf, 'GenCT-Cog_beh_overlay.fig')
saveas(gcf, 'GenCT-Cog_beh_overlay.png')
close all

%% FIGURE 3 - SpecCT-Cog LV1

% PREP
load pls1_vertex_SpecCT-Cog\civetsurf_pls_SpecCT-Cog.mat
pls1 = PLS;
avsurf = data1.avsurf;
mask = data1.mask;
load pls2_vertex_SpecCT-Cog\civetsurf_pls_SpecCT-Cog.mat
pls2 = PLS;
clearvars -except pls1 pls2 avsurf mask behdesc

% SURFSTAT OVERLAY 
U1 = pls1.result.boot_result.compare_u(:,1); % comp1
U2 = pls2.result.boot_result.compare_u(:,1); % comp1

%pos
U1pos = U1;
U1pos(U1pos<0) = 0;
U2pos = U2;
U2pos(U2pos<0) = 0;
surfstat_overlay(U1pos,U2pos,mask,avsurf)
SurfStatColLim([4,7])
cmap = [0.8 0.8 0.8; autumn(3)];
SurfStatColormap(cmap(1:4,:))
saveas(gcf, 'SpecCT-Cog_overlay_pos.fig')
saveas(gcf, 'SpecCT-Cog_overlay_pos.png')
close all

%neg
U1neg = U1;
U1neg(U1neg>0) = 0;
U2neg = U2;
U2neg(U2neg>0) = 0;
surfstat_overlay(U1neg,U2neg,mask,avsurf)
SurfStatColLim([1,4])
cmap = [0.8 0.8 0.8; cool(3)];
SurfStatColormap(cmap)
saveas(gcf, 'SpecCT-Cog_overlay_neg.fig')
saveas(gcf, 'SpecCT-Cog_overlay_neg.png')
close all

% BEH OVERLAY
pls_overlay (pls1, 1, pls2, 1, behdesc, [1,1])
saveas(gcf, 'SpecCT-Cog_beh_overlay.fig')
saveas(gcf, 'SpecCT-Cog_beh_overlay.png')
close all

%% Save glimfiles w/ component scores

clear
clc
%sample1
load pls1_vertex_GenCT-Cog\civetsurf_pls_GenCT-Cog.mat
d1 = d1p.glimfile;
d1.genCTcog1 = PLS.result.usc(:,1);
clearvars -except PLS d1

load('pls1_vertex_SpecCT-Cog\civetsurf_pls_SpecCT-Cog.mat', 'PLS')
d1.specCTcog1 = PLS.result.usc(:,1);
d1.specCTcog2 = PLS.result.usc(:,2)*-1; % flipped

writetable(d1, 'data1_glimfile_pls.csv')

clear
clc
load pls2_vertex_GenCT-Cog\civetsurf_pls_GenCT-Cog.mat
d2 = d2p.glimfile;
d2.genCTcog1 = PLS.result.usc(:,1);
clearvars -except PLS d2

load('pls2_vertex_SpecCT-Cog\civetsurf_pls_SpecCT-Cog.mat', 'PLS')
d2.specCTcog1 = PLS.result.usc(:,1);
d2.specCTcog2 = PLS.result.usc(:,2)*-1; % flipped
d2.specCTcog3 = PLS.result.usc(:,3)*-1; % flipped

writetable(d2, 'data2_glimfile_pls.csv')