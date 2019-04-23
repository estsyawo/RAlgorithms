# create full correlation matrix from lower triangular matrix

nvars <- 12

#-------------------------------------------------------->

# the following is going to be read columnwise even though it looks rowwise
lowerFacCorMatSym <- c(1.000,-0.732,-0.357,0.752,-0.568,-0.402,0.644,-0.487,-0.319,0.618,-0.501,-0.269,
                       1.000,0.476,-0.534,0.689,0.460,-0.455,0.622,0.406,-0.377,0.543,0.287,
                       1.000,-0.352,0.423,0.616,-0.349,0.417,0.612,-0.221,0.309,0.563,
                       1.000,-0.714,-0.474,0.649,-0.456,-0.362,0.602,-0.408,-0.247,
                       1.000,0.624,-0.516,0.738,0.570,-0.424,0.550,0.350,
                       1.000,-0.404,0.530,0.680,-0.272,0.425,0.625,
                       1.000,-0.647,-0.464,0.543,-0.530,-0.352,
                       1.000,0.663,-0.338,0.694,0.433,
                       1.000,-0.347,0.490,0.569,
                       1.000,-0.555,-0.326,
                       1.000,0.575,
                       1.000)
# create the full correlation matrix 
X <- diag(nvars)
X[lower.tri(X, diag=TRUE)] <- lowerFacCorMatSym
X <- X + t(X) - diag(diag(X)) 
longFacCorMatFull <- X
longFacCorMatFull