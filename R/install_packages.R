inst_packages <- installed.packages()
for (package in c("tidyverse", "sf", "raster", "rmapshaper", "spdep", 
                  "spatstat",  "mapview", "fpc", "GWmodel", "NbClust", 
                  "reshape2", "mapedit", "osmdata", "patchwork", "terra", 
                  "exactextractr", "fasterize")) {
    test <- package %in% inst_packages[,1]
    if (!test) {
        install.packages(package)
    }
}
rm(package, inst_packages, test)