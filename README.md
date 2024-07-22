# ProcessMicroplate
Matlab Live Script file for processing bacterial growth data from an automated microtiter plate reader. Users require Matlab and a working knowledge of how to open and run Matlab LiveScript files. The MPmodel3.m Matlab class file must be in the Matlab path.  

The files include:
1. EX240513.csv - example raw data (microtiter platereader output file) for processing
2. EX240513_C1C2.csv - example out file for columns 1 and 2
3. EX240513_C3C4.csv - example out file for columns 3 and 4
4. EX240513_C5C6.csv - example out file for columns 5 and 6
5. EX240513_C7C8.csv - example out file for columns 7 and 8
6. EX240513_C9C10.csv - example out file for columns 9 and 10
7. MPmodel3.m - Matlab class with processing functions
8. Processplate.mlx - Matlab LiveScript file for processing microplate data
9. LICENSE - CC0 license

Matlab is required to run the LiveScript program, and MPmodel3.m must be in the the Matlab path. 

ProcessMicroplate
 A Matlab Live-Script Program 

Version: 3.6, 2024.06.28

Author: Fred Breidt, USDA/ARS Research Microbiologist, Food Science and Market Quality and Handling Research Unit, 322 Schaub Hall, Box 7624, NC State University, Raleigh, NC 27695-7624
Email: fred.breidt@usda.gov

Introduction
While software exists for determining microbial growth kinetics from 96 well microtiter plate optical density data, the accuracy of current methods has been questioned. Variation in derived parameters can occur due to improper background subtraction and failing to consider the impact of initial optical density on the derived parameters (Atolia et al., 2020; Begot et al., 1996). 

Background

The MPmodel3 class has the flexibility to allow the user to address these deficiencies and observe the influence initial conditions has on the kinetic parameter output, including:
Growth rate (GR)
Lag time (LAG)
Initial or minimum optical density (MIN)
Maximum optical density (MAX)
Note that GR, LAG, MIN, and MAX are output table variable names for each of the growth kinetics parameters (see below).

The MPmodel3 class has functions to obtain growth parameters using a sliding window method (Breidt, 1994; Atolia et al., 2020) that is not depending on modeling a sigmoidal shaped curve. The advantage of this approach is that growth parameters can be accurately obtained from oddly shaped growth curves that don’t look like the traditions sigmoidal curve. An additional feature of the MPmodel3 class is that the output parameters can be used to generate a predicted Gompertz curve using a published algorithm (Zewitering et al. 1990).      

Table 1. Data input format. This file is easy to generate and required only minor modification (using Microsoft Excel) from the original microtiter plate reader output .csv file. The file is saved as a .csv file and then opened using this script (see Step 1 of the LiveScript file). 

For an example, see file EX240513.csv

Table 2. Output data. Data is output to .csv tables with the text name as entered (see Step 2 of the LiveScript file). The output vairable and corresponding standard deviation for duplicates (i.e. for columns 1 and 2 of the example plate data).

For an example, see file EX240513_C1C2.csv 
 
The output data files may be opened directly in MS Excel (as .csv file) for further processing. Graphs may be copied to Microsoft Word, PowerPoint or other software using the MS Windows clipboard.   

References: 

Atolia E, Cesar S, Arjes HA, Rajendram M, Shi H, Knapp BD, Khare S, Aranda-Díaz A, Lenski RE, Huang KC. 2020. Environmental and physiological factors affecting high-throughput measurements of bacterial growth. mBio. 11(5):e01378-20. https://doi.org/10.1128/mBio.01378-201

Begot C, Desnier I, Daudin JD, Labadie JC, Lebert A. 1996. Recommendations for calculating growth parameters by optical density measurements. J Microbiol Meth 25(3):225-232. https://doi.org/10.1016/0167-7012(95)00090-9

Breidt F, Romick TL, Fleming HP. 1994. A rapid method for the determination of bacterial growth kinetics. J Rapid Methods Autom Microbiol 3(1):59-68. https://doi.org/10.1111/j.1745-4581.1994.tb00308.x

Zwietering MH, Jongenburger I, Rombouts FM, van't Tiet K. 1990. Modeling the bacterial growth curve. Appl Env Microbiol 56(6):1875-1881. https://doi.org/10.1128/aem.56.6.1875-1881.1990
