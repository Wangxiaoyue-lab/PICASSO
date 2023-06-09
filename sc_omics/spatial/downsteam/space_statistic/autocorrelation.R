# moran's I
library(spdep)
library(spatialreg)

listw <- nb2listw(your_neighbours, style = "W")
moran.test(your_data, listw)


# Geary's C
geary.test(your_data, listw)





# Getis and Ord's G
getisord(your_data, listw)
