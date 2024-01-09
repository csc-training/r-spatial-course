# Repository for CSC 'Spatial data analysis with R' course 

The aim of this course is to familiarize participants with spatial analysis with R. 
All materials by Marko Kallio (Aalto University).
Last instructor led course was on 10-12-5.2023: https://ssl.eventilla.com/event/rspatial_23

## Contents
### Day 1, vector basics
* D1S1 Handling and plotting vector data in R
* D1S2 Handling and plotting cont.
* D1S3 Spatial operations (intersection, clipping, conversions etc)
* D1S4 Spatial operations cont.
 
### Day 2, vector data analysis and visualization
* D2S1 Spatial analysis of vector data (clustering, density surfaces, autocorrelation)
* D2S2 Spatial analysis cont.
* D2S3 Visualizing spatial data

### Day 3, raster basics
* D2S1 Raster basics with R
* D3S2 Raster data manipulation
* D3S3 Map algebra
* D3S4 Spatial modelling with raster data

### Self-study recommendations
* Vector data exercises (Days 1  and 2) should be done in the given order, with possibility to skip to visualization before finishing some of the sessions before it.
* Raster data exercises (Day 3) can be done as own module, without doing the vector data exercises.
* Use RStudio to open the material, each session has its own R notebook.
* Each session ends with exercises, sample solutions can be found from separate file for each day.
* To open first session open in RStudio the folder `r-spatial-course` and then to specific day.
* The "Visual" mode of RStudio notebooks seems to sometimes have problems, use then the "Source" mode.
* There are two main options for self-study, use CSC Noppe service or with local RStudio.

#### Using CSC Noppe

CSC Noppe has self-study RStudio with ready package installations and course materials. It is available for Finnish academic users with HAKA or VIRTU account.

* Open [CSC Noppe](https://notebooks.rahtiapp.fi/) as exercise environment.
* Start `Spatial data analysis with R` public application, which has a copy of this repository.
* CSC Noppe has rather limited memory, so after each session restart R from Session menu.
* The public notebook does not save any files at session end, so export any files you want before finishing.

#### Using local RStudio

* Install RStudio
* Install [required R packages](install_packages.sh)
* Get materials. Clone [this Github repository](https://github.com/csc-training/r-spatial-course). In RStudio: `File -> New project -> Version control -> Git`
  * Repository URL: `https://github.com/csc-training/r-spatial-course.git`
  * Project directory name: `r-spatial-course`

## Prerequisites
* Basics of geoinformatics and geostatistics
* Basic use of R, no earlier experience with R spatial packages is needed. For self-study [Data analysis with R course material](https://github.com/csc-training/da-with-r-remote) can be used.

## Additional links
* [Other external learning materials](https://docs.csc.fi/apps/r-env-for-gis/#references)
