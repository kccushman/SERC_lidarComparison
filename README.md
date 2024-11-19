# Case study: Multi-scale forest ecological inference from lidar platforms

This repository includes data and code accompanying the case study included in the manuscript "Multi-scale forest ecological inference from lidar platform" (submitted to _Ecology Letters_). We compiled lidar datasets from multiple platforms in a common area to:

1.	Demonstrate how differences in sensor characteristics influence density and resolution of lidar data.
2.	Provide open-source, co-located datasets for users to further inspect differences in lidar data.
3.	Provide example code to perform basic lidar analysis.

This case study is meant to allow readers to get hands-on experience with real-world data from different platforms. This case study is not meant to be a rigorous comparison of derived ecological metrics among all sensors; such comparisons can be found throughout other publications referenced throughout the main manuscript. Code includes basic functions in R commonly used to visualize and manipulate lidar data accessible with a normal laptop computer; more sophisticated algorithms for advanced users are also referenced throughout the main manuscript.

Terrestrial laser scanning (TLS), mobile laser scanning (MLS), UAS laser scanning (ULS), airborne laser scanning (ALS), and spaceborne laser scanning (SLS) data were collected within the Smithsonian Environmental Research Center (SERC) forest dynamics plot in Maryland, USA. TLS, MLS, and ALS data were collected within 1 month of the 2021 growing season; ULS data were collected in November 2020 (“leaf-off” data) and July 2022 (“leaf-on” data). 

To facilitate visualization and comparison of discrete return lidar clouds, we chose an 80 x 5 m transect through the center of the study area. This transect is in the center of the MLS trajectory, the dataset with the smallest spatial extent (Fig. S1). For each instrument, we calculated voxelized point density, gridded canopy height (1 m grid), total and voxelized effective plant area index (ePAI), and canopy rugosity (spatial variation in vertical variation of voxelized ePAI; Hardiman et al. 2011) within the transect. We estimated ePAI using the MacArthur-Horn approach, which estimates ePAI from the proportion of point measurements (lasers) that encounter vegetation or pass through a voxel (MacArthur and Horn, 1969). We use the term “ePAI” instead of “leaf area index” because we did not attempt to adjust ePAI estimates based on sensor or vegetation structural properties (clumping, leaf angle, wood) to illustrate differences in sampling properties among platforms. For ULS, ePAI and canopy rugosity were only calculated using leaf-on data.

## Data folder
- "GEDI" folder: GEDI data from the initial ~ 4 year deployment over the SERC NEON study area. NOTE: some other GEDI folders need to be downloaded from the following Google Drive link https://drive.google.com/drive/folders/11y4RBRfPJSHqWiekDEiomRX5UlqNJeF0?usp=sharing . Keep organization in "GEDI01_B", "GEDI02_A", "GEDI02_B" folders within this GEDI folder.
- "ha4_data" ALS, ULS, and MLS data from ha 4 in the SERC plot. NOTE: TLS data need to be downloaded from the following Google Drive link:
 https://drive.google.com/drive/folders/1zTAWz94pnOWlbFfPlgONtNq5IFj58dgW?usp=sharing and unzipped in a "TLS" folder.
- "transect" ALS, ULS, MLS, and TLS files are included for the 80 x 5 m transect within the SERC plot ha 4. The subsetting is performed in Code/SubsetClouds.R
- "trunk" manually subsetted cross-section of a single tree trunk in TLS, MLS, and ULS (leaf off) data

## Code folder
- GEDI_openData.R : convert GEDI .h5 files (as originally downloaded) to summarized .csv files of metrics needed for analysis
- GEDI_analyzeData.R : analyze summarized GEDI .csv files to produce results and figures for the manuscript
- SubsetClouds.R : subset 80 x 5 m transect from larger TLS, MLS, ULS, and ALS point clouds from ha 4 of the SERC plot
- TransectStats.R : analyze 80 x 5 m transect data from TLS, MLS, ULS, and ALS point clouds to produce results and figures for the manuscript

## Results folder
Contains data-based figures are results produced from the included Data and Code, as numbered in the submitted manuscript

## Other resources
In addition to the R code included here, we encourage readers to vizualize point clouds in 3D using CloudeCompare, a free and open source software: https://cloudcompare-org.danielgm.net/ 

Manuscript authors: K.C. Cushman (MLS data), Angelica M. Almeyda Zambrano (ULS data), Eben N. Broadbent (ULS data), Kamil Král (TLS data), Martin Krůček (TLS data), Sean M. McMahon (SERC plot PI). Data analysis and visualization by KCC, contact: cushmankc@ornl.gov .
