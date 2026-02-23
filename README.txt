Data and code for the study "Localised negative feedback shapes genome-wide patterning of meiotic DNA breaks"

Instructions for the simulator:
-Ensure R version 4.4 or later is installed. Rstudio installation is also recommended.
-Load and run the DSBint_function script. This will register the simulator as a memory object.
-Call in the simulator by using DSBint(). Hotspot tables need to be supplied to the WTmap (map with interference) and tel1Dmap (map with low/no interference, used as the template for the simulation). These can either be R memory objects or paths to hotspot table files (see data folder AVERAGE_HOTSPOT_TABLES for examples). Other parameters have defaults that may be changed. See the DSBint_function script for more details on parameters.
-Once finished, the simulator will consolidate all runs into a single output file. The column "sim" contains the simulated hits, in NormHpChr. "WT" denotes the number of hits in the WT hotspot table, while "tel1D" is the number of hits in the template hotspot table. To compare input and simulator data, use the columns "SimRatio" and "RealRatio".

The simulator has been tested on Mac and Linux. It is dependent on the R packages data.table, plyr, e1071, tictoc, bio3d, doParallel, foreach, doFuture, dplyr, and bspec, all of which it will attempt to install upon its first usage if not installed prior to running.

Scripts:
R SCRIPTS: Averaging_FullMap_tables_V1: This script averages individual FullMap replicates into a combined FullMap where the sum of HpM equals 1 million.
Calculating background reads_V1: This script measures the percentage of signal registered within the 50 largest genes—regions of presumed Spo11 inactivity—on the S. cerevisiae genome as an estimate of the background noise.
DSBint_function: This script is the DSB simulator. The simulator generates virtual chromosomes of user-specified length and number, binning coordinates by user-specified bin widths. Each bin is assigned a relative probability based on the input DSB hotspot map. For more details on potential parameters, see the help for the function.
Hotspot_identification_V1: This script identifies position and length of hotspots on single or multiple Spo11-DSB libraries. The total HpM signal is smoothed with a 201 Hann window. A cutoff of 0.1 HpM is then applied to remove the background noise. Hotspots are defined setting a minimum length of 25 bp and a minimum number of reads of 25. Hotspots separated by < 200 bp are merged and considered as a single hotspot. Hotspots are defined in each library separately and then combined to produce a single hotspot template that defines the position of every hotspot identified on the libraries.
Hotspot_Smooth_ratios_V1:This script calculates and represents the hotspots' fold changes between two Libraries (NormHpM or NormHpChr ratio).
Hotspot_table_V1: This script calculates the HpM signal included within each hotspot. Detailed description of the term reference list included in “Hotspot Table Definitions_01”.

Data:
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
AVERAGE_FULLMAPs:
FullMap.Cer3H4L2_MJ906_sae2rec114-8A_Average_RA32_RA76.txt
FullMap.Cer3H4L2_MJ1439_sae2xrs2-KF_Average_RA53_RA54.txt
FullMap.Cer3H4L2_MJ1438_sae2xrs2-11_Average_RA51_RA52.txt
FullMap.Cer3H4L2_LLR36_sae2ndt80rec1148D_Average_RA45_RA46.txt
FullMap.Cer3H4L2_LLR40_sae2ndt80rec1148Dtel1_AverageRA47_RA48.txt 
FullMap.Cer3H4L2_MJ1099_sae2Dnd80Dtel1-kd_Average_1A_ccLLR16.txt
FullMap.Cer3H4L2_MJ907_sae2Drec114-8D_Average_RA33_ccLLR15.txt
FullMap.Cer3H4L2_MJ962_sae2Dndt80D_AverageD1D2TC10TC17.txt
FullMap.Cer3H4L2_MJ965_sae2Dndt80Dtel1D_AverageD1D2TC5TC10TC17.txt
FullMap.Cer3H4L2_sae2Dndt80D_SPO11-GBDtel1D_AverageccLLR13-21.txt
FullMap.Cer3H4L2_sae2Dndt80D_SPO11-GBD_Average_ccLLR14-22.txt
FullMap.Cer3H4L2_VG402_sae2Dtel1D_Average_1AB234.txt
FullMap.Cer3H4L2_MJ315_sae2D_Average_12A357BCtxt.txt

Pseudocode version of simulator script:

FUNCTION DSBInt(parameters)

1. Set Working Directory
   IF running in RStudio:
       set working directory to script location
   ELSE IF script path available:
       set working directory to script location

2. Load Required Packages
   FOR each package in package_list:
       IF package not installed:
           install package
       load package

3. Define Helper: load_hotspot(x)
   IF x is dataframe:
       RETURN x
   ELSE IF x is valid file path:
       read file into dataframe
       RETURN dataframe
   ELSE:
       STOP with error

4. Load Hotspot Data
   A.hotspot ← load_hotspot(tel1Dmap)
   B.hotspot ← load_hotspot(WTmap)

   IF "NormHpMChr" exists AND "NormHpChr" missing:
       copy NormHpMChr → NormHpChr

   Keep columns: Chr, Midpoint, NormHpChr

5. Compute Global Parameters
   E1 ← max(RepN)
   RepTN ← round(G / RepN)
   maxsim ← max(RepTN)

6. Precompute Exponential Windows (if Winmethod == "Exponential")
   FOR each window width w in Windows:
       compute exponential window of width w
       store in exp_windows dictionary

7. Build Output Directory Name
   Generate timestamp string
   IF fail_on == TRUE:
       format fail_label with failrate
   ELSE:
       fail_label ← empty

   fname ← formatted string including:
       ActiveChromatids
       G
       RepN range
       Winmethod
       Windows range
       chrom
       fail_label
       timestamp
       script_name

   Create directories:
       output_dir/sim_parts/(optional expname)/fname
       output_dir/plots/fname
       output_dir/sim_plots/fname

8. Subset Chromosome Data
   A1 ← subset A.hotspot for chromosome chrom
   B1 ← subset B.hotspot for chromosome chrom

   Convert Midpoint to bin index (round(Midpoint/res))

   Aggregate duplicate bins by summing

   L ← round(ChrSizes[chrom] / res)

   Initialize arrays:
       A[1..L] ← 0
       B[1..L] ← 0

   Fill A using A1 NormHpChr values
   Fill B using B1 NormHpChr values

9. Determine Chromosome DSB Count
   E_raw ← scaled DSB count based on chromosome length and ActiveChromatids
   E ← round(E_raw)
   IF E < 1:
       E ← 1

   Repseq ← floor((E/E1) * RepN)
   Initialize AllResult ← NULL

10. Parallel Loop Over Window Sizes
   PARALLEL FOR each window k in Windows:

       FOR each transStrength tr:

           W ← round(k*1000/res)
           constrain W between 1 and 2*L

           IF Winmethod == "Hann":
               C ← 1 - hanning.window(W)
           ELSE IF Winmethod == "Exponential":
               C ← 1 - exp_windows[k]
           ELSE:
               C ← 1 - tukeywindow(W, winr)

           Initialize:
               H[1..E] ← zero arrays length L
               F[1..E] ← 0
               cellcount ← empty table

11. Simulate Cells
           FOR i = 1 to maxsim:

               DSBN ← 0
               DSBS ← 0
               DSBF ← 0
               RepIndex ← 1

               FOR each chromatid d in 1..ActiveChromatids:
                   D[d][1..L] ← 1

               FOR j = 1 to E:

                   IF i <= RepTN[RepIndex]:

                       Dact ← random chromatid

                       weighted_probs ← A * D[Dact]
                       IF sum(weighted_probs) == 0:
                           CONTINUE

                       pos ← sample bin using weighted_probs

                       is_failed ← FALSE
                       IF fail_on == TRUE:
                           IF random() < failrate:
                               is_failed ← TRUE

                       DSBN++

                       IF is_failed:
                           DSBF++
                       ELSE:
                           DSBS++
                           H[j][pos]++

                       F[j]++

12. Build Interference Window Around pos
                       left ← max(1, pos - floor(W/2))
                       right ← min(L, pos + floor(W/2))

                       Extract appropriate portion of C
                       Adjust length if necessary
                       C1[left:right] ← window slice

                       scale_factor ← failStrength IF failed ELSE 1

13. Apply Trans Interference
                       IF tr > 0:
                           FOR each chromatid ≠ Dact:
                               D[chromatid] ← D[chromatid] *
                                   (((C1 * scale_factor) + ((1/tr)-1)) * tr)
                               enforce minimum threshold

14. Apply Cis Interference
                       D[Dact] ← D[Dact] *
                           (((C1 * scale_factor) + ((1/cisStrength)-1)) * cisStrength)
                       enforce minimum threshold

15. Snapshot at Replication Checkpoints
                       IF j in Repseq:
                           success_frac ← DSBS / DSBN
                           success_count ← round(success_frac * RepN[RepIndex])
                           failed_count ← RepN[RepIndex] - success_count

                           append to cellcount:
                               DSBs = DSBN
                               GDSBs = RepN[RepIndex]
                               Cells = i
                               success_frac
                               success_count
                               failed_count

                           RepIndex++

16. Post-Simulation Aggregation
           H1[1] ← H[1]
           FOR h = 2 to E:
               H1[h] ← H1[h-1] + H[h]

           FOR each h:
               IF sum(H1[h]) > 0:
                   normalize H1[h] to per million
               ELSE:
                   set H1[h] to zero vector

17. Generate Output Tables
           FOR each h in Repseq:

               Create Result table:
                   Pos
                   tel1D = A
                   WT = B
                   sim = H1[h]

               Compute:
                   SimRatio = log2(sim/tel1D)
                   RealRatio = log2(WT/tel1D)

               Remove NA and infinite rows

               Smooth:
                   SimR
                   RealR
                   DevR

               rmsd ← sqrt(mean((SimRatio - RealRatio)^2))
               rmsdA ← sqrt(mean(DevR^2))

               cellsub ← subset cellcount where GDSBs matches h

               Add metadata columns:
                   method
                   window
                   N.attempted
                   successful
                   failed
                   success_frac
                   transstrength

               Write Result to CSV in sdir

18. Merge Simulation Parts
   files ← list files in sdir
   AllResult ← concatenate all CSV files
   delete individual part files

   Write AllResult to final CSV

19. Summarize Results
   FOR each method:
       FOR each window:
           FOR each transstrength:
               compute mean successful DSBs
               compute mean failed DSBs
               append to simsum table

20. Print "Sim generation done"

END FUNCTION
