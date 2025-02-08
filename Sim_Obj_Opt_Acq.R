rm(list=ls())
require(FITSio)
require(fields)
require(ggplot2)
require(plotrix)

out_Ref <- "~/Desktop/Pol-Gal/Simulation/Reference/"
out_Beams <- "~/Desktop/Pol-Gal/Simulation/Beams/"
out_Pol <- "~/Desktop/Pol-Gal/Simulation/Polarization/"
out_Res <- "~/Desktop/Pol-Gal/Simulation/Residuals/"

#Size of the object to simulate
N <- 20
n <- N / 2
K <- N * N
k <- n * n

#Generating matrices to hold simulated object
Ref_Obj <- array(NaN, dim = c(N, N, 7),
                 dimnames = list(NULL, NULL, c("I", "Ip", "Iu", "P", 
                                               "Q", "U", "X")))
E_Field_Obj <- array(NaN, dim = c(N, N, 3), 
                     dimnames = list(NULL, NULL, c("E", "Ex", "Ey")))
Quad1_Obj <- array(NaN, dim = c(n, n, 6),
                   dimnames = list(NULL, NULL, c("I", "Ip", "Iu", 
                                                 "Q", "U", "X")))
Quad2_Obj <- array(NaN, dim = c(n, n, 6),
                   dimnames = list(NULL, NULL, c("I", "Ip", "Iu", 
                                                 "Q", "U", "X")))
Quad3_Obj <- array(NaN, dim = c(n, n, 6),
                   dimnames = list(NULL, NULL, c("I", "Ip", "Iu", 
                                                 "Q", "U", "X")))
Quad4_Obj <- array(NaN, dim = c(n, n, 6),
                   dimnames = list(NULL, NULL, c("I", "Ip", "Iu", 
                                                 "Q", "U", "X")))

PR_raw_Obj <- array(NaN, dim = c(N, N, 4),
                    dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))
PR_IF_Obj <- array(NaN, dim = c(N, N, 4),
                   dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))
PR_PF_Obj <- array(NaN, dim = c(N, N, 4),
                   dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))
PR_IPF_Obj <- array(NaN, dim = c(N, N,4),
                    dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))

B_raw_Obj <- array(NaN, dim = c(N, N,4),
                   dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))
B_IF_Obj <- array(NaN, dim = c(N, N, 4),
                  dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))
B_PF_Obj <- array(NaN, dim = c(N, N, 4),
                  dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))
B_IPF_Obj <- array(NaN, dim = c(N, N, 4),
                   dimnames = list(NULL, NULL, c("P", "Q", "U", "X")))

#Filling Quadrant matrices
# Quad1_Obj[, , "I"] <- 1
# Quad2_Obj[, , "I"] <- .75
# Quad3_Obj[, , "I"] <- .5
# Quad4_Obj[, , "I"] <- .25
Quad1_Obj[, , "I"] <- 1000
Quad2_Obj[, , "I"] <- 1000
Quad3_Obj[, , "I"] <- 1000
Quad4_Obj[, , "I"] <- 1000

P <- 1
# for(i in 1:n){
#   for(j in 1:n){
#     Quad1_Obj[i, j, "Ip"] <- (1 - (i * j - 1) / n^2) * Quad1_Obj[i, j, "I"]
#     Quad2_Obj[i, j, "Ip"] <- (1 - (i * j - 1) / n^2) * Quad2_Obj[i, j, "I"]
#     Quad3_Obj[i, j, "Ip"] <- (1 - (i * j - 1) / n^2) * Quad3_Obj[i, j, "I"]
#     Quad4_Obj[i, j, "Ip"] <- (1 - (i * j - 1) / n^2) * Quad4_Obj[i, j, "I"]
#   }
# }


Quad1_Obj[, , "Ip"] <- P * Quad1_Obj[, , "I"]
Quad2_Obj[, , "Ip"] <- P * Quad2_Obj[, , "I"]
Quad3_Obj[, , "Ip"] <- P * Quad3_Obj[, , "I"]
Quad4_Obj[, , "Ip"] <- P * Quad4_Obj[, , "I"]

Quad1_Obj[, , "Iu"] <- (1 - P) * Quad1_Obj[, , "I"]
Quad2_Obj[, , "Iu"] <- (1 - P) * Quad2_Obj[, , "I"]
Quad3_Obj[, , "Iu"] <- (1 - P) * Quad3_Obj[, , "I"]
Quad4_Obj[, , "Iu"] <- (1 - P) * Quad4_Obj[, , "I"]

# for(i in 1:n){
#   for(j in 1:n){
#     Quad1_Obj[i, j, "Ip"] <- Ipol_r * Quad1_Obj[i, j, "I"]
#     Quad2_Obj[i, j, "Ip"] <- Ipol_r * Quad2_Obj[i, j, "I"]
#     Quad3_Obj[i, j, "Ip"] <- Ipol_r * Quad3_Obj[i, j, "I"]
#     Quad4_Obj[i, j, "Ip"] <- Ipol_r * Quad4_Obj[i, j, "I"]
#   }
# }

Xs <- seq(0, 1/2, length.out = K) #cospi(x) <=> cos(pi*x)

# Quad1_Obj[, , "X"] <- matrix(Xs[1:k], nrow = n, byrow = FALSE)
# Quad2_Obj[, , "X"] <- matrix(Xs[(k + 1):(2 * k)], nrow = n, byrow = FALSE)
# Quad3_Obj[, , "X"] <- matrix(Xs[(2 * k + 1):(3 * k)], nrow = n, byrow = FALSE)
# Quad4_Obj[, , "X"] <- matrix(Xs[(3 * k + 1):(4 * k)], nrow = n, byrow = FALSE)
Quad1_Obj[, , "X"] <- matrix(rep(0, k), nrow = n, byrow = FALSE)
Quad2_Obj[, , "X"] <- matrix(rep(1/8, k), nrow = n, byrow = FALSE)
Quad3_Obj[, , "X"] <- matrix(rep(1/4, k), nrow = n, byrow = FALSE)
Quad4_Obj[, , "X"] <- matrix(rep(3/8, k), nrow = n, byrow = FALSE)

Quad1_Obj[, , "Q"] <- Quad1_Obj[, , "Ip"] * cospi(2 * Quad1_Obj[, , "X"])
Quad2_Obj[, , "Q"] <- Quad2_Obj[, , "Ip"] * cospi(2 * Quad2_Obj[, , "X"])
Quad3_Obj[, , "Q"] <- Quad3_Obj[, , "Ip"] * cospi(2 * Quad3_Obj[, , "X"])
Quad4_Obj[, , "Q"] <- Quad4_Obj[, , "Ip"] * cospi(2 * Quad4_Obj[, , "X"])

Quad1_Obj[, , "U"] <- Quad1_Obj[, , "Ip"] * sinpi(2 * Quad1_Obj[, , "X"])
Quad2_Obj[, , "U"] <- Quad2_Obj[, , "Ip"] * sinpi(2 * Quad2_Obj[, , "X"])
Quad3_Obj[, , "U"] <- Quad3_Obj[, , "Ip"] * sinpi(2 * Quad3_Obj[, , "X"])
Quad4_Obj[, , "U"] <- Quad4_Obj[, , "Ip"] * sinpi(2 * Quad4_Obj[, , "X"])

#Carrying Quadrant matrices into tensor
Quads_Obj <- array(NaN, dim = c(n, n, 6, 4))
Quads_Obj[, , , 1] <-Quad1_Obj
Quads_Obj[, , , 2] <-Quad2_Obj
Quads_Obj[, , , 3] <-Quad3_Obj
Quads_Obj[, , , 4] <-Quad4_Obj
rm(Quad1_Obj, Quad2_Obj, Quad3_Obj, Quad4_Obj)

#Filling Unit matrix with quadrant matrices
for(i in 1:4){
  switch(i,
         c(j <- 1, k <- n, l <- 1, m <- n),
         c(j <- n + 1, k <- 2 * n, l <- 1, m <- n),
         c(j <- n + 1, k <- 2 * n, l <- n + 1, m <- 2 * n),
         c(j <- 1, k <- n, l <-n + 1, m <- 2 * n))

  Ref_Obj[j:k, l:m, -4] <- Quads_Obj[, , , i]
}
rm(Quads_Obj)

Ref_Obj[, , "P"] <- Ref_Obj[, , "Ip"] / Ref_Obj[, , "I"]

#Filling E field matrix
E_Field_Obj[, , "E"] <- sqrt(2 * Ref_Obj[, , "I"])
unpol_E <- sqrt(2 * Ref_Obj[, , "Iu"])
unpol_Ex <- unpol_E / sqrt(2)
unpol_Ey <- unpol_Ex

pol_E <- sqrt(2 * Ref_Obj[, , "Ip"])
pol_Ex <- pol_E * cospi(Ref_Obj[, , "X"])
pol_Ey <- pol_E * sinpi(Ref_Obj[, , "X"])

E_Field_Obj[, , "Ex"] <- unpol_Ex + pol_Ex
E_Field_Obj[, , "Ey"] <- unpol_Ey + pol_Ey
  
#Create FITS files for Ref I, P, X
writeFITSim(Ref_Obj[, , "I"], file = paste0(out_Ref, "Ref_I.fits"))
writeFITSim(Ref_Obj[, , "P"], file = paste0(out_Ref, "Ref_P.fits"))
writeFITSim(Ref_Obj[, , "X"] * 180, file = paste0(out_Ref, "Ref_X.fits"))
writeFITSim(Ref_Obj[, , "Q"], file = paste0(out_Ref, "Ref_Q.fits"))
writeFITSim(Ref_Obj[, , "U"], file = paste0(out_Ref, "Ref_U.fits"))

pdf(paste0(out_Ref, "Ref-arrows.pdf"), width = 37, height = 32)
par(mfrow = c(1, 1))
par(mar = c(15, 15, 15, 40))
par(cex.main = 6)
par(cex.axis = 6)
par(cex.lab = 6)

temp_x <- 0:(N - 1)
temp_x <- rep(temp_x, N) / (N - 1)
temp_y <- 0:(N - 1)
temp_y <- rep(temp_y, each = N) / (N - 1)

temp_P <- Ref_Obj[, , "P"]
temp_X <- as.vector(Ref_Obj[, , "X"]) * 180
temp_M <- as.vector(temp_P)
temp_scale <- 8 * max(temp_P) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Simulation Reference",
  legend.width = 3, legend.mar = 0
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")
dev.off()

print("Reference Object has been simulated.")

#Generating I flat field
O_flat_gen <- rnorm(K, mean = 0, sd = 0.3)
O_flat_gen <- O_flat_gen - 1.01 * min(O_flat_gen)
O_flat_gen <- O_flat_gen / max(O_flat_gen)
O_I_flat <- matrix(O_flat_gen, nrow = N)

E_flat_gen <- rnorm(K, mean = 0, sd = 0.3)
E_flat_gen <- E_flat_gen - 1.01 * min(E_flat_gen)
E_flat_gen <- E_flat_gen / max(E_flat_gen)
E_I_flat <- matrix(E_flat_gen, nrow = N)

#Generating Pol flat
O_P_flat <- 0.2
E_P_flat <- 0.7

#Signal goes through Half-Wave Plate
pos_ang <- c(0, 1/8, 1/4, 3/8)
E_Beam <- array(NaN, dim = c(N, N, 4),
                dimnames = list(NULL, NULL, c("0", "22", "45", "67")))
O_Beam <- array(NaN, dim = c(N, N, 4),
                dimnames = list(NULL, NULL, c("0", "22", "45", "67")))

for(l in 1:4){
  O_Beam[, , l] <- .5 * Ref_Obj[, , "Iu"] +
    Ref_Obj[, , "Ip"] * cospi(2 * pos_ang[l] - Ref_Obj[, , "X"])^2

  E_Beam[, , l] <- .5 * Ref_Obj[, , "Iu"] +
    Ref_Obj[, , "Ip"] * sinpi(2 * pos_ang[l] - Ref_Obj[, , "X"])^2
}

Ref_Obj[, , "Q"] <- Ref_Obj[, , "Q"] / Ref_Obj[, , "I"]
Ref_Obj[, , "U"] <- Ref_Obj[, , "U"] / Ref_Obj[, , "I"]
Ref_Obj[, , "X"] <- Ref_Obj[, , "X"] * 180
Ref_Obj <- Ref_Obj[, , -(1:3)]

writeFITSim(E_Beam, file = paste0(out_Beams, "E_Beams_HWP.fits"))
writeFITSim(O_Beam, file = paste0(out_Beams, "O_Beams_HWP.fits"))

print("The signal has passed through the HWP.")

#Signal goes through Wollaston Prism
E_Beam <- E_Beam * E_P_flat
O_Beam <- O_Beam * O_P_flat

writeFITSim(E_Beam, file = paste0(out_Beams, "E_Beams_WP.fits"))
writeFITSim(O_Beam, file = paste0(out_Beams, "O_Beams_WP.fits"))

print("The signal has passed through the WP.")

#Signal reaches detector chip
for(b in 1:4){
  E_Beam[, , b] <- E_Beam[, , b] * E_I_flat
  O_Beam[, , b] <- O_Beam[, , b] * O_I_flat
}

writeFITSim(E_Beam, file = paste0(out_Beams, "E_Beams_CCD.fits"))
writeFITSim(O_Beam, file = paste0(out_Beams, "O_Beams_CCD.fits"))

print("The signal has reached the CCD.")

#Measured, uncalibrated Beams are processed
#using both Patat&Romaniello and Bagnulo's method
beam_comb_PR_raw <- (O_Beam - E_Beam) / (O_Beam + E_Beam)
beam_comb_B_raw <- O_Beam / E_Beam

PR_raw_Obj[, , "Q"] <- .5 * (beam_comb_PR_raw[,,"0"] - beam_comb_PR_raw[,,"45"])
PR_raw_Obj[, ,"U"] <- .5 * (beam_comb_PR_raw[,,"22"] - beam_comb_PR_raw[,,"67"])
PR_raw_Obj[, , "P"] <- sqrt(PR_raw_Obj[, , "Q"]^2 + PR_raw_Obj[, , "U"]^2)
temp_X <- .5 * atan(PR_raw_Obj[, , "U"] / PR_raw_Obj[, , "Q"]) * 180 /pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
PR_raw_Obj[, , "X"] <- temp_X
rm(beam_comb_PR_raw, temp_ind, temp_X)

temp_a <- beam_comb_B_raw[, , "0"] / beam_comb_B_raw[, , "45"]
temp_b <- beam_comb_B_raw[, , "22"] / beam_comb_B_raw[, , "67"]
B_raw_Obj[, , "Q"] <- (sqrt(temp_a) - 1) / (sqrt(temp_a) + 1)
B_raw_Obj[, , "U"] <- (sqrt(temp_b) - 1) / (sqrt(temp_b) + 1)
B_raw_Obj[, , "P"] <- sqrt(B_raw_Obj[, , "Q"]^2 + B_raw_Obj[, , "U"]^2)
temp_X <- .5 * atan(B_raw_Obj[, , "U"] / B_raw_Obj[, , "Q"]) * 180 / pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
B_raw_Obj[, , "X"] <- temp_X
rm(beam_comb_B_raw, temp_a, temp_b, temp_ind, temp_X)

writeFITSim(PR_raw_Obj[, , "Q"], file = paste0(out_Pol, "PR_raw_Q.fits"))
writeFITSim(PR_raw_Obj[, , "U"], file = paste0(out_Pol, "PR_raw_U.fits"))
writeFITSim(PR_raw_Obj[, , "P"], file = paste0(out_Pol, "PR_raw_P.fits"))
writeFITSim(PR_raw_Obj[, , "X"], file = paste0(out_Pol, "PR_raw_X.fits"))
writeFITSim(B_raw_Obj[, , "Q"], file = paste0(out_Pol, "B_raw_Q.fits"))
writeFITSim(B_raw_Obj[, , "U"], file = paste0(out_Pol, "B_raw_U.fits"))
writeFITSim(B_raw_Obj[, , "P"], file = paste0(out_Pol, "B_raw_P.fits"))
writeFITSim(B_raw_Obj[, , "X"], file = paste0(out_Pol, "B_raw_X.fits"))

print("Polarization has been estimated without correcting for I or Pol.")

#Measured, I flatted beams are processed
#using both Patat&Romaniello and Bagnulo's method
temp_O <- O_Beam
temp_E <- E_Beam
for(b in 1:4){
  temp_O[, , b] <- temp_O[, , b] / O_I_flat
  temp_E[, , b] <- temp_E[, , b] / E_I_flat
}
beam_comb_PR_IF <- (temp_O - temp_E) / (temp_O + temp_E)
beam_comb_B_IF <- temp_O / temp_E
rm(temp_E, temp_O)

PR_IF_Obj[, , "Q"] <- .5 * (beam_comb_PR_IF[,,"0"] - beam_comb_PR_IF[,,"45"])
PR_IF_Obj[, , "U"] <- .5 * (beam_comb_PR_IF[,,"22"] - beam_comb_PR_IF[,,"67"])
PR_IF_Obj[, , "P"] <- sqrt(PR_IF_Obj[, , "Q"]^2 + PR_IF_Obj[, , "U"]^2)
temp_X <- .5 * atan(PR_IF_Obj[, , "U"] / PR_IF_Obj[, , "Q"]) * 180 / pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
PR_IF_Obj[, , "X"] <- temp_X
rm(beam_comb_PR_IF, temp_ind, temp_X)

temp_a <- beam_comb_B_IF[, , "0"] / beam_comb_B_IF[, , "45"]
temp_b <- beam_comb_B_IF[, , "22"] / beam_comb_B_IF[, , "67"]
B_IF_Obj[, , "Q"] <- (sqrt(temp_a) - 1) / (sqrt(temp_a) + 1)
B_IF_Obj[, , "U"] <- (sqrt(temp_b) - 1) / (sqrt(temp_b) + 1)
B_IF_Obj[, , "P"] <- sqrt(B_IF_Obj[, , "Q"]^2 + B_IF_Obj[, , "U"]^2)
temp_X <- .5 * atan(B_IF_Obj[, , "U"] / B_IF_Obj[, , "Q"]) * 180 / pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
B_IF_Obj[, , "X"] <- temp_X
rm(beam_comb_B_IF, temp_a, temp_b, temp_ind, temp_X)

writeFITSim(PR_IF_Obj[, , "Q"], file = paste0(out_Pol, "PR_IF_Q.fits"))
writeFITSim(PR_IF_Obj[, , "U"], file = paste0(out_Pol, "PR_IF_U.fits"))
writeFITSim(PR_IF_Obj[, , "P"], file = paste0(out_Pol, "PR_IF_P.fits"))
writeFITSim(PR_IF_Obj[, , "X"], file = paste0(out_Pol, "PR_IF_X.fits"))
writeFITSim(B_IF_Obj[, , "Q"], file = paste0(out_Pol, "B_IF_Q.fits"))
writeFITSim(B_IF_Obj[, , "U"], file = paste0(out_Pol, "B_IF_U.fits"))
writeFITSim(B_IF_Obj[, , "P"], file = paste0(out_Pol, "B_IF_P.fits"))
writeFITSim(B_IF_Obj[, , "X"], file = paste0(out_Pol, "B_IF_X.fits"))

print("Polarization has been estimated correcting for I but not for Pol.")

#Measured, P flatted beams are processed
#using both Patat&Romaniello and Bagnulo's method
temp_O <- O_Beam / O_P_flat
temp_E <- E_Beam / E_P_flat
beam_comb_PR_PF <- (temp_O - temp_E) / (temp_O + temp_E)
beam_comb_B_PF <- temp_O / temp_E
rm(temp_E, temp_O)

PR_PF_Obj[, , "Q"] <- .5 * (beam_comb_PR_PF[,,"0"] - beam_comb_PR_PF[,,"45"])
PR_PF_Obj[, , "U"] <- .5 * (beam_comb_PR_PF[,,"22"] - beam_comb_PR_PF[,,"67"])
PR_PF_Obj[, , "P"] <- sqrt(PR_PF_Obj[, , "Q"]^2 + PR_PF_Obj[, , "U"]^2)
temp_X <- .5 * atan(PR_PF_Obj[, , "U"] / PR_PF_Obj[, , "Q"]) * 180 / pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
PR_PF_Obj[, , "X"] <- temp_X
rm(beam_comb_PR_PF, temp_ind, temp_X)

temp_a <- beam_comb_B_PF[, , "0"] / beam_comb_B_PF[, , "45"]
temp_b <- beam_comb_B_PF[, , "22"] / beam_comb_B_PF[, , "67"]
B_PF_Obj[, , "Q"] <- (sqrt(temp_a) - 1) / (sqrt(temp_a) + 1)
B_PF_Obj[, , "U"] <- (sqrt(temp_b) - 1) / (sqrt(temp_b) + 1)
B_PF_Obj[, , "P"] <- sqrt(B_PF_Obj[, , "Q"]^2 + B_PF_Obj[, , "U"]^2)
temp_X <- .5 * atan(B_PF_Obj[, , "U"] / B_PF_Obj[, , "Q"]) * 180 / pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
B_PF_Obj[, , "X"] <- temp_X
rm(beam_comb_B_PF, temp_a, temp_b, temp_ind, temp_X)

writeFITSim(PR_PF_Obj[, , "Q"], file = paste0(out_Pol, "PR_PF_Q.fits"))
writeFITSim(PR_PF_Obj[, , "U"], file = paste0(out_Pol, "PR_PF_U.fits"))
writeFITSim(PR_PF_Obj[, , "P"], file = paste0(out_Pol, "PR_PF_P.fits"))
writeFITSim(PR_PF_Obj[, , "X"], file = paste0(out_Pol, "PR_PF_X.fits"))
writeFITSim(B_PF_Obj[, , "Q"], file = paste0(out_Pol, "B_PF_Q.fits"))
writeFITSim(B_PF_Obj[, , "U"], file = paste0(out_Pol, "B_PF_U.fits"))
writeFITSim(B_PF_Obj[, , "P"], file = paste0(out_Pol, "B_PF_P.fits"))
writeFITSim(B_PF_Obj[, , "X"], file = paste0(out_Pol, "B_PF_X.fits"))

print("Polarization has been estimated correcting for Pol but not for I.")

#Measured, I and P flatted beams are processed
#using both Patat&Romaniello and Bagnulo's method
temp_O <- O_Beam / O_P_flat
temp_E <- E_Beam / E_P_flat
for(b in 1:4){
  temp_O[, , b] <- temp_O[, , b] / O_I_flat
  temp_E[, , b] <- temp_E[, , b] / E_I_flat
}
beam_comb_PR_IPF <- (temp_O - temp_E) / (temp_O + temp_E)
beam_comb_B_IPF <- temp_O / temp_E
rm(temp_E, temp_O)

PR_IPF_Obj[,, "Q"] <- .5 * (beam_comb_PR_IPF[,, "0"] - beam_comb_PR_IPF[,,"45"])
PR_IPF_Obj[,, "U"] <- .5 * (beam_comb_PR_IPF[,,"22"] - beam_comb_PR_IPF[,,"67"])
PR_IPF_Obj[,, "P"] <- sqrt(PR_IPF_Obj[, , "Q"]^2 + PR_IPF_Obj[, , "U"]^2)
temp_X <- .5 * atan(PR_IPF_Obj[, , "U"] / PR_IPF_Obj[, , "Q"]) * 180 / pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
PR_IPF_Obj[,, "X"] <- temp_X
rm(beam_comb_PR_IPF, temp_ind, temp_X)

temp_a <- beam_comb_B_IPF[, , "0"] / beam_comb_B_IPF[, , "45"]
temp_b <- beam_comb_B_IPF[, , "22"] / beam_comb_B_IPF[, , "67"]
B_IPF_Obj[, , "Q"] <- (sqrt(temp_a) - 1) / (sqrt(temp_a) + 1)
B_IPF_Obj[, , "U"] <- (sqrt(temp_b) - 1) / (sqrt(temp_b) + 1)
B_IPF_Obj[, , "P"] <- sqrt(B_IPF_Obj[, , "Q"]^2 + B_IPF_Obj[, , "U"]^2)
temp_X <- .5 * atan(B_IPF_Obj[, , "U"] / B_IPF_Obj[, , "Q"]) * 180 / pi
temp_ind <- which(temp_X < 0, arr.ind = TRUE)
temp_X[temp_ind] <- temp_X[temp_ind] + 90
B_IPF_Obj[, , "X"] <- temp_X
rm(beam_comb_B_IPF, temp_a, temp_b, temp_ind, temp_X)

writeFITSim(PR_IPF_Obj[, , "Q"], file = paste0(out_Pol, "PR_IPF_Q.fits"))
writeFITSim(PR_IPF_Obj[, , "U"], file = paste0(out_Pol, "PR_IPF_U.fits"))
writeFITSim(PR_IPF_Obj[, , "P"], file = paste0(out_Pol, "PR_IPF_P.fits"))
writeFITSim(PR_IPF_Obj[, , "X"], file = paste0(out_Pol, "PR_IPF_X.fits"))
writeFITSim(B_IPF_Obj[, , "Q"], file = paste0(out_Pol, "B_IPF_Q.fits"))
writeFITSim(B_IPF_Obj[, , "U"], file = paste0(out_Pol, "B_IPF_U.fits"))
writeFITSim(B_IPF_Obj[, , "P"], file = paste0(out_Pol, "B_IPF_P.fits"))
writeFITSim(B_IPF_Obj[, , "X"], file = paste0(out_Pol, "B_IPF_X.fits"))

print("Polarization has been estimated correcting for both I and Pol.")

#Polarization Arrow Plots for each case
pdf(paste0(out_Pol, "PR-arrows.pdf"), width = 75, height = 65)
par(mfrow = c(2, 2))
par(mar = c(15, 15, 15, 20))
par(cex.main = 6)
par(cex.axis = 6)
par(cex.lab = 6)

xs <- 0:(N - 1)
xs <- rep(xs, N) / (N - 1)
ys <- 0:(N - 1)
ys <- rep(ys, each = N) / (N - 1)

temp_P <- PR_raw_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(PR_raw_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(PR_raw_Obj[, , "X"])[-temp_inds]
}

temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Patat&Romaniello, RAW",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")

temp_P <- PR_IF_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(PR_IF_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(PR_IF_Obj[, , "X"])[-temp_inds]
}

temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Patat&Romaniello, IF",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")

temp_P <- PR_PF_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(PR_PF_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(PR_PF_Obj[, , "X"])[-temp_inds]
}

temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Patat&Romaniello, PF",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")

temp_P <- PR_IPF_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(PR_IPF_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(PR_IPF_Obj[, , "X"])[-temp_inds]
}

temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Patat&Romaniello, IPF",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")
dev.off()

pdf(paste0(out_Pol, "B-arrows.pdf"), width = 75, height = 75)
par(mfrow = c(2, 2))
par(mar = c(15, 15, 15, 20))
par(cex.main = 6)
par(cex.axis = 6)
par(cex.lab = 6)

temp_P <- B_raw_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(B_raw_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(B_raw_Obj[, , "X"])[-temp_inds]
}

temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Bagnulo, RAW",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")

temp_P <- B_IF_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(B_IF_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(B_IF_Obj[, , "X"])[-temp_inds]
}

temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Bagnulo, IF",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")

temp_P <- B_PF_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(B_PF_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(B_PF_Obj[, , "X"])[-temp_inds]
}

temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Bagnulo, PF",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")

temp_P <- B_IPF_Obj[, , "P"]
temp_data <- as.vector(temp_P)
temp_inds <- which(is.na(temp_data), arr.ind = TRUE)
rm(temp_data)

temp_x <- xs
temp_y <- ys
temp_M <- as.vector(temp_P)
temp_X <- as.vector(B_IPF_Obj[, , "X"])

if(length(temp_inds) != 0){
  temp_x <- xs[-temp_inds]
  temp_y <- ys[-temp_inds]
  temp_M <- as.vector(temp_P)[-temp_inds]
  temp_X <- as.vector(B_IPF_Obj[, , "X"])[-temp_inds]
}
temp_scale <- 8 * max(temp_P, na.rm = TRUE) / N

image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Polarization Map, Bagnulo, IPF",
  legend.width = 6, legend.mar = 1
)
vectorField(temp_X, temp_M, temp_x, temp_y, scale = temp_scale,
            headspan = 0.01, vecspec = "deg", col = "red")
dev.off()
rm(temp_M, temp_P, temp_X, temp_x, temp_y)

print("Polarization files have been created.")

#Residuals for each polarization estimate are calculated
temp_Res <- abs((PR_raw_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(PR_raw_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
PR_raw_Res <- temp_Res

temp_Res <- abs((B_raw_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(B_raw_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
B_raw_Res <- temp_Res

temp_Res <- abs((PR_IF_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(PR_IF_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
PR_IF_Res <- temp_Res

temp_Res <- abs((B_IF_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(B_IF_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
B_IF_Res <- temp_Res

temp_Res <- abs((PR_PF_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(PR_PF_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
PR_PF_Res <- temp_Res

temp_Res <- abs((B_PF_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(B_PF_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
B_PF_Res <- temp_Res

temp_Res <- abs((PR_IPF_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(PR_IPF_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
PR_IPF_Res <- temp_Res

temp_Res <- abs((B_IPF_Obj - Ref_Obj) / Ref_Obj) * 100
check_inds <- which(is.na(temp_Res) | is.infinite(temp_Res), arr.ind = TRUE)
check_length <- dim(check_inds)[1]
if(!is.null(check_length)){
  for(i in 1:check_length){
    c <- check_inds[i, ]
    if(Ref_Obj[c[1], c[2], c[3]] == 0){
      temp_Res[c[1], c[2], c[3]] <- abs(B_IPF_Obj[c[1], c[2], c[3]]) * 100
    }else{
      temp_Res[c[1], c[2], c[3]] <- NA
    }
  }
}
B_IPF_Res <- temp_Res

rm(B_IF_Obj, B_IPF_Obj, B_PF_Obj, B_raw_Obj,
   PR_IF_Obj, PR_IPF_Obj, PR_PF_Obj, PR_raw_Obj, temp_Res)

writeFITSim(PR_raw_Res[, , "Q"], file = paste0(out_Res, "PR_raw_Q_Res.fits"))
writeFITSim(PR_raw_Res[, , "U"], file = paste0(out_Res, "PR_raw_U_Res.fits"))
writeFITSim(PR_raw_Res[, , "P"], file = paste0(out_Res, "PR_raw_P_Res.fits"))
writeFITSim(PR_raw_Res[, , "X"], file = paste0(out_Res, "PR_raw_X_Res.fits"))
writeFITSim(B_raw_Res[, , "Q"], file = paste0(out_Res, "B_raw_Q_Res.fits"))
writeFITSim(B_raw_Res[, , "U"], file = paste0(out_Res, "B_raw_U_Res.fits"))
writeFITSim(B_raw_Res[, , "P"], file = paste0(out_Res, "B_raw_P_Res.fits"))
writeFITSim(B_raw_Res[, , "X"], file = paste0(out_Res, "B_raw_X_Res.fits"))
writeFITSim(PR_IF_Res[, , "Q"], file = paste0(out_Res, "PR_IF_Q_Res.fits"))
writeFITSim(PR_IF_Res[, , "U"], file = paste0(out_Res, "PR_IF_U_Res.fits"))
writeFITSim(PR_IF_Res[, , "P"], file = paste0(out_Res, "PR_IF_P_Res.fits"))
writeFITSim(PR_IF_Res[, , "X"], file = paste0(out_Res, "PR_IF_X_Res.fits"))
writeFITSim(B_IF_Res[, , "Q"], file = paste0(out_Res, "B_IF_Q_Res.fits"))
writeFITSim(B_IF_Res[, , "U"], file = paste0(out_Res, "B_IF_U_Res.fits"))
writeFITSim(B_IF_Res[, , "P"], file = paste0(out_Res, "B_IF_P_Res.fits"))
writeFITSim(B_IF_Res[, , "X"], file = paste0(out_Res, "B_IF_X_Res.fits"))
writeFITSim(PR_PF_Res[, , "Q"], file = paste0(out_Res, "PR_PF_Q_Res.fits"))
writeFITSim(PR_PF_Res[, , "U"], file = paste0(out_Res, "PR_PF_U_Res.fits"))
writeFITSim(PR_PF_Res[, , "P"], file = paste0(out_Res, "PR_PF_P_Res.fits"))
writeFITSim(PR_PF_Res[, , "X"], file = paste0(out_Res, "PR_PF_X_Res.fits"))
writeFITSim(B_PF_Res[, , "Q"], file = paste0(out_Res, "B_PF_Q_Res.fits"))
writeFITSim(B_PF_Res[, , "U"], file = paste0(out_Res, "B_PF_U_Res.fits"))
writeFITSim(B_PF_Res[, , "P"], file = paste0(out_Res, "B_PF_P_Res.fits"))
writeFITSim(B_PF_Res[, , "X"], file = paste0(out_Res, "B_PF_X_Res.fits"))
writeFITSim(PR_IPF_Res[, , "Q"], file = paste0(out_Res, "PR_IPF_Q_Res.fits"))
writeFITSim(PR_IPF_Res[, , "U"], file = paste0(out_Res, "PR_IPF_U_Res.fits"))
writeFITSim(PR_IPF_Res[, , "P"], file = paste0(out_Res, "PR_IPF_P_Res.fits"))
writeFITSim(PR_IPF_Res[, , "X"], file = paste0(out_Res, "PR_IPF_X_Res.fits"))
writeFITSim(B_IPF_Res[, , "Q"], file = paste0(out_Res, "B_IPF_Q_Res.fits"))
writeFITSim(B_IPF_Res[, , "U"], file = paste0(out_Res, "B_IPF_U_Res.fits"))
writeFITSim(B_IPF_Res[, , "P"], file = paste0(out_Res, "B_IPF_P_Res.fits"))
writeFITSim(B_IPF_Res[, , "X"], file = paste0(out_Res, "B_IPF_X_Res.fits"))

pdf(paste0(out_Res, "PR-Res.pdf"), width = 85, height = 90)
par(mfrow = c(2, 2))
par(mar = c(15, 15, 15, 20))
par(cex.main = 6)
par(cex.axis = 6)
par(cex.lab = 6)

temp_P <- PR_raw_Res[, , "P"]
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Patat&Romaniello, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_P <- PR_IF_Res[, , "P"]
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Patat&Romaniello, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_P <- PR_PF_Res[, , "P"]
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Patat&Romaniello, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_P <- PR_IPF_Res[, , "P"]
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Patat&Romaniello, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- PR_raw_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Patat&Romaniello, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- PR_IF_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Patat&Romaniello, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- PR_PF_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Patat&Romaniello, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- PR_IPF_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Patat&Romaniello, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- PR_raw_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Patat&Romaniello, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- PR_IF_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Patat&Romaniello, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- PR_PF_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Patat&Romaniello, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- PR_IPF_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Patat&Romaniello, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- PR_raw_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Patat&Romaniello, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- PR_IF_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Patat&Romaniello, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- PR_PF_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Patat&Romaniello, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- PR_IPF_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Patat&Romaniello, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

dev.off()

pdf(paste0(out_Res, "B-Res.pdf"), width = 85, height = 90)
par(mfrow = c(2, 2))
par(mar = c(15, 15, 15, 20))
par(cex.main = 6)
par(cex.axis = 6)
par(cex.lab = 6)

temp_P <- B_raw_Res[, , "P"]
temp_x <- 1:N
temp_y <- 1:N
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Bagnulo, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_P <- B_IF_Res[, , "P"]
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Bagnulo, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_P <- B_PF_Res[, , "P"]
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Bagnulo, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_P <- B_IPF_Res[, , "P"]
image.plot(
  temp_P,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for P, Bagnulo, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_P, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- B_raw_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Bagnulo, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- B_IF_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Bagnulo, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- B_PF_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Bagnulo, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_X <- B_IPF_Res[, , "X"]
image.plot(
  temp_X,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for X, Bagnulo, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_X, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- B_raw_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Bagnulo, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- B_IF_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Bagnulo, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- B_PF_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Bagnulo, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_Q <- B_IPF_Res[, , "Q"]
image.plot(
  temp_Q,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for Q, Bagnulo, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_Q, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- B_raw_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Bagnulo, RAW",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- B_IF_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Bagnulo, IF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- B_PF_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Bagnulo, PF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")

temp_U <- B_IPF_Res[, , "U"]
image.plot(
  temp_U,
  xlab='Rows', ylab='Cols', col = gray.colors(128, 0, 1, 2.6),
  xaxt="n", yaxt="n", main = "Residuals for U, Bagnulo, IPF",
  legend.width = 6, legend.mar = 1
)
mtext(paste0("Mean residual: ", mean(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 8, cex = 5, col = "red")
mtext(paste0("Median residual: ", median(temp_U, na.rm = TRUE), "%"), 
      side = 1, line = 14, cex = 5, col = "Blue")
dev.off()

print("Residuals files have been created.")