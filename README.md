Data and code for the study "Localised negative feedback shapes genome-wide patterning of meiotic DNA breaks"

Instructions for the simulator:
-Ensure R version 4.4 or later is installed. Rstudio installation is also recommended.
-Load and run the DSBint_function script. This will register the simulator as a memory object.
-Call in the simulator by using DSBint(). Hotspot tables need to be supplied to the WTmap (map with interference) and tel1Dmap (map with low/no interference, used as the template for the simulation). These can either be R memory objects or paths to hotspot table files (see data folder AVERAGE_HOTSPOT_TABLES for examples). Other parameters have defaults that may be changed. See the DSBint_function script for more details on parameters.
-Once finished, the simulator will consolidate all runs into a single output file. The column "sim" contains the simulated hits, in NormHpChr. "WT" denotes the number of hits in the WT hotspot table, while "tel1D" is the number of hits in the template hotspot table. To compare input and simulator data, use the columns "SimRatio" and "RealRatio".

The simulator has been tested on Mac and Linux. It is dependent on the R packages data.table, plyr, e1071, tictoc, bio3d, doParallel, foreach, doFuture, dplyr, and bspec, all of which it will attempt to install upon its first usage if not installed prior to running.

HS_TEMPLATE: Averages Hotspot template used in this studied. Hotspots identified in rDNA region were excluded.
MGBD_Hotspot_template_4019HS 
T_0.193_Hotspot_template_3473HS 
AVERAGE_HOTSPOT_TABLES: Averages of the biological replicates, plus hotspot template used to measure hotspot-specific signals in each library: 
T_0.193_Hotspot.Table.sae2D.txt 
T_0.193_Hotspot.Table.sae2Dtel1D.txt
T_0.193_Hotspot.Table.sae2Dndt80D.txt 
T_0.193_Hotspot.Table.sae2Dndt80Dtel1D.txt
T_0.193_Hotspot.Table.sae2Dndt80Dtel1kd.txt
T_0.193_Hotspot.Table.sae2Dxrs2-11.txt
T_0.193_Hotspot.Table.sae2Dxrs2-KF.txt
T_0.193_Hotspot.Table.sae2Drec114.8A.txt
T_0.193_Hotspot.Table.sae2Drec114.8D.txt
T_0.193_Hotspot.Table.sae2Dndt80Drec1148D.txt
T_0.193_Hotspot.Table.sae2Dndt80Drec1148Dtel1D.txt
MGBD_Hotspot.Table.SPO11-GBD.txt
MGBD_Hotspot.Table.SPO11-GBDtel1D.txt

Scripts:
R SCRIPTS: Averaging_FullMap_tables_V1: This script averages individual FullMap replicates into a combined FullMap where the sum of HpM equals 1 million.
Calculating background reads_V1: This script measures the percentage of signal registered within the 50 largest genes—regions of presumed Spo11 inactivity—on the S. cerevisiae genome as an estimate of the background noise.
DSBint_function: This script is the DSB simulator. The simulator generates virtual chromosomes of user-specified length and number, binning coordinates by user-specified bin widths. Each bin is assigned a relative probability based on the input DSB hotspot map. For more details on potential parameters, see the help for the function.
Hotspot_identification_V1: This script identifies position and length of hotspots on single or multiple Spo11-DSB libraries. The total HpM signal is smoothed with a 201 Hann window. A cutoff of 0.1 HpM is then applied to remove the background noise. Hotspots are defined setting a minimum length of 25 bp and a minimum number of reads of 25. Hotspots separated by < 200 bp are merged and considered as a single hotspot. Hotspots are defined in each library separately and then combined to produce a single hotspot template that defines the position of every hotspot identified on the libraries.
Hotspot_Smooth_ratios_V1:This script calculates and represents the hotspots' fold changes between two Libraries (NormHpM or NormHpChr ratio).
Hotspot_table_V1: This script calculates the HpM signal included within each hotspot. Detailed description of the term reference list included in “Hotspot Table Definitions_01”.
