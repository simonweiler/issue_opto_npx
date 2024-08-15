This is the data structure and code for the manuscript: 
Overcoming off-target optical stimulation-evoked cortical activity in the mouse brain in vivo
Simon Weiler, Mateo Velez-Fort  and Troy W. Margrie

Data_structure: Data_00-25-Jul-2024.mat 

animal_name

exp_name

genotype (1 = wild type) 

session_ids (1-12 sessions) 
%1) blue darkness (1-5 laser intensities)
%2) orange darkness (1-5 laser intensities)
%3) red darkness (1-5 laser intensities)

%4) red ambient light 20 lux (1-5 laser intensities)
%5) red ambient light 40 lux (1-5 laser intensities)
%6) red ambient light 80 lux (1-5 laser intensities)

%7) blue ambient light 20 lux (1-5 laser intensities)
%8) blue ambient light 40 lux (1-5 laser intensities)
%9) blue ambient light 80 lux (1-5 laser intensities)

%10) orange ambient light 20 lux (1-5 laser intensities)
%11) orange ambient light 40 lux (1-5 laser intensities)
%12) orange ambient light 80 lux (1-5 laser intensities)

sessiontype (105-116): random identifier for a given session

include (1 = include; 0= exclude)

spike_type (1=FS; 2=RS)

spike_type_dur (duration of spike to idnetify FS and RS)

session_aligned
	raster: spike time data for all 50 stimulations for SINGLE unit activity (10 trials per laser intensity, 5 laser intensities) 
	reaster_mua: spike time data for all 50 stimulations for MULTI unit activity (10 trials per laser intensity, 5 laser intensities) 
	sdf_cortex: Spike time density for each SINGLE unit for each laser intensity (cells x duration x nr intensity)
	sdf_mua: Spike time density for each MULTI unit for each laser intensity (cells x duration x nr intensity)
	param_cortex: Modulaiton index and delta firing rate for each cell (cells x MI/dFR x nr intensity)
	zeta_cortex: ZETA test result across intensities
	ttest_cortex: ttest result across intensities
	zeta_cortex_noise: ZETA test result across intensities for CONTROL
	ttest_cortex_noise: ttest result across intensities for CONTROL
	param_mua: Modulaiton index and delta firing rate for each MULTI UNIT (cells x MI/dFR x nr intensity)
	zeta_mua: ZETA test result across intensities for each MULTI UNIT
	ttest_mua: ttest result across intensities for each MULTI UNIT








