* Encoding: UTF-8.
CROSSTABS
  /TABLES=group BY study
  /FORMAT=AVALUE TABLES
  /CELLS=COUNT COLUMN BPROP 
  /COUNT ROUND CELL.

* Group differences - Combined samples.
CROSSTABS
  /TABLES= sex handedness BY group
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN BPROP 
  /COUNT ROUND CELL.
  
  T-TEST GROUPS=Group(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=age education iq z_ps z_att z_wm z_vm z_vism z_ef z_gci mean_thickness20mm_2_1_1_harmonized
  /ES DISPLAY(TRUE)
  /CRITERIA=CI(.95).

*Check battery.
 UNIANOVA z_ps BY group WITH battery
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(group) WITH(battery=MEAN) 
  /PRINT ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=battery group.

UNIANOVA z_att BY group WITH battery
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(group) WITH(battery=MEAN) 
  /PRINT ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=battery group.

UNIANOVA z_wm BY group WITH battery
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(group) WITH(battery=MEAN) 
  /PRINT ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=battery group.

UNIANOVA z_vm BY group WITH battery
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(group) WITH(battery=MEAN) 
  /PRINT ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=battery group.

UNIANOVA z_vism BY group WITH battery
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(group) WITH(battery=MEAN) 
  /PRINT ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=battery group.

UNIANOVA z_ef BY group WITH battery
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(group) WITH(battery=MEAN) 
  /PRINT ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=battery group.

UNIANOVA z_gci BY group WITH battery
  /METHOD=SSTYPE(3)
  /INTERCEPT=INCLUDE
  /EMMEANS=TABLES(group) WITH(battery=MEAN) 
  /PRINT ETASQ
  /CRITERIA=ALPHA(.05)
  /DESIGN=battery group.

*Sample 1.
temporary.
select if study = 1.
CROSSTABS
  /TABLES=sex handedness BY group
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN BPROP 
  /COUNT ROUND CELL.

temporary.
select if study = 1.
T-TEST GROUPS=Group(1 2)
/MISSING=ANALYSIS
/VARIABLES=age education iq z_ps z_att z_wm z_vm z_vism z_ef z_gci mean_thickness20mm_2_1_1_harmonized
/ES DISPLAY(TRUE)
/CRITERIA=CI(.95).

temporary.
select if study = 1.
DESCRIPTIVES VARIABLES= saps sans cpz_adherence dur_illness
  /STATISTICS=MEAN STDDEV MIN MAX.

VALUE LABELS dx 295.10 'schizophrenia-disorganized' 295.20 'schizophrenia-catatonic' 295.30 'schizophrenia-paranoid' 295.40 'schizophreniform' 295.60 'schizophrenia - residual type' 295.70 'schizoaffective' 
    295.90 'schizophrenia, undifferentiated type' 296.00 'bipolar I disorder, single manic episode, unspecified' 296.04 'bipolar I with psychotic features' 296.34 'major depressive disorder, recurrent, severe with psychotic features'
    296.40 'bipolar disorder, manic, unspecified' 296.44 'bipolar disorder, manic, with psychotic features' 296.50 'bipolar disorder, depressed, unspecified' 296.64 'bipolar 1, mixed, severe, with psychotic features' 
    296.89 'bipolar II' 297.10 'delusional disorder' 298.90 'psychotic disorder, NOS' 301.20 'schizoid personality disorder' 292.10 'substance-induced', 296.80 'bipolar, unspecified', 298.8 'brief psychotic disorder'.
    
temporary.
select if study = 1 & group = 2.
FREQUENCIES VARIABLES= dx
  /ORDER=ANALYSIS.

*Sample 2.
temporary.
select if study = 2.
CROSSTABS
  /TABLES= sex handedness BY group
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN BPROP 
  /COUNT ROUND CELL.

temporary.
select if study = 2.
T-TEST GROUPS=Group(1 2)
/MISSING=ANALYSIS
/VARIABLES=age education iq z_ps z_att z_wm z_vm z_vism z_ef z_gci mean_thickness20mm_2_1_1_harmonized
/ES DISPLAY(TRUE)
/CRITERIA=CI(.95).

temporary.
select if study = 2.
DESCRIPTIVES VARIABLES= saps sans cpz_adherence dur_illness
  /STATISTICS=MEAN STDDEV MIN MAX.

temporary.
select if study = 2 & Group = 2.
FREQUENCIES VARIABLES = dx
  /ORDER=ANALYSIS.

*Combined samples - Patient.
temporary.
select if Group = 2.
CROSSTABS
  /TABLES= sex handedness BY study
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN BPROP 
  /COUNT ROUND CELL.

temporary.
select if Group = 2.
T-TEST GROUPS=Study(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=age education iq z_ps z_att z_wm z_vm z_vism z_ef z_gci mean_thickness20mm_2_1_1_harmonized saps sans cpz_adherence dur_illness
  /CRITERIA=CI(.95).

*Combined samples - Control.
temporary.
select if Group = 1.
CROSSTABS
  /TABLES= sex handedness BY study
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ 
  /CELLS=COUNT COLUMN BPROP 
  /COUNT ROUND CELL.

temporary.
select if Group = 1.
T-TEST GROUPS=Study(1 2)
  /MISSING=ANALYSIS
  /VARIABLES=age education iq z_ps z_att z_wm z_vm z_vism z_ef z_gci mean_thickness20mm_2_1_1_harmonized
  /CRITERIA=CI(.95).
