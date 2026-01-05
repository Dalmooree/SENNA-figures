suppressPackageStartupMessages(
  suppressWarnings({
    library(SENNA)
    library(Seurat)
    library(dplyr)
    library(foreach)
    library(doParallel)
}))

# path
dpath <- paste0(getwd(), "/data/Thymus9/outs")

## Utility functions----
progsimulation <- function(signal_intensity,
                           n = 1,
                           poipar,
                           senna) {
  cprog <- senna
  ref <- cprog@Coord[["Spatial"]]
  ref[["t"]] <- ref[["t"]] - min(ref[["t"]])
  dist <- max(ref[["distance"]]) - ref[["distance"]]
  dist <- (dist / max(dist))^2
  t <- ref[["t"]] / max(ref[["t"]])
  
  
  grad_index <- c(signal_intensity * poipar[1:500],
                  - signal_intensity * poipar[501:1000],
                  rep(0, 1e3))
  params <- outer(t * dist, grad_index) + 
    matrix(rep(poipar, nrow(ref)), byrow = TRUE, nrow = nrow(ref))

  params[params < 0] <- 0
  
  save_dir <- file.path(getwd(), "data", "sim.generated", "sim.prog")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  set.seed(seed = n)
  prog_count <- t(
    matrix(
      rpois(length(params), lambda = params), 
      nrow = nrow(params),
      dimnames = list(rownames(ref),
                      colnames(cprog@Gene[["Spatial"]]))))
  
  write.csv(prog_count,
            file.path(save_dir, paste0("sim_prog_SI", round(signal_intensity, 2), ".csv")),
            row.names = TRUE,
            quote = TRUE)

  prog_perm <- t(
    matrix(
      apply(prog_count, 1, sample),
      ncol = nrow(prog_count),
      dimnames = list(colnames(prog_count),
                      rownames(prog_count))))

  write.csv(prog_perm,
            file.path(save_dir, paste0("perm_prog_SI", round(signal_intensity, 2), ".csv")),
            row.names = TRUE,
            quote = TRUE)
  
  # Ground truth
  gs <- rownames(prog_count)[1:1000]
  gc <- rownames(prog_count)[1001:2000]
  
  # 기존 코드와의 호환성을 위해 wide 포맷으로 저장
  ground_truth <- cbind(gs, gc)
  colnames(ground_truth) <- c("gs", "gc")
  write.csv(ground_truth,
            file.path(save_dir, paste0("GT_SI", round(signal_intensity, 2), ".csv")), 
            row.names = FALSE, 
            quote = TRUE)
}


regiosimulation <- function(signal_intensity,
                            n = 1,
                            poipar,
                            senna) {
  cregio <- senna
  ref <- cregio@Coord[["Spatial"]]
  csd <- (ref[["distance"]] - min(ref[["distance"]])) / 
    (max(ref[["distance"]]) - min(ref[["distance"]]))
  csd[csd == 0] <- 1e-10
  
  grad_index <- c(signal_intensity * poipar[1:500],
                  - signal_intensity * poipar[501:1000],
                  rep(0, 1e3))
  params <- outer(csd, grad_index) + 
    matrix(rep(poipar, nrow(ref)), byrow = TRUE, nrow = nrow(ref)) 

  params[params < 0] <- 0

  save_dir <- file.path(getwd(), "data", "sim.generated", "sim.regio")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  set.seed(seed = n)
  regio_count <- t(
    matrix(
      rpois(length(params), lambda = params), 
      nrow = nrow(params),
      dimnames = list(rownames(ref),
                      colnames(cregio@Gene[["Spatial"]]))))

  write.csv(regio_count,
            file.path(save_dir, paste0("sim_regio_SI", round(signal_intensity, 2), ".csv")),
            row.names = TRUE,
            quote = TRUE)
  
  regio_perm <- t(
    matrix(
      apply(regio_count, 1, sample),
      ncol = nrow(regio_count),
      dimnames = list(colnames(regio_count),
                      rownames(regio_count))))
  write.csv(regio_perm,
            file.path(save_dir, paste0("perm_regio_SI", round(signal_intensity, 2), ".csv")),
            row.names = TRUE,
            quote = TRUE)
  
  
  # Ground truth
  gs <- rownames(regio_count)[1:1000]
  gc <- rownames(regio_count)[1001:2000]

  # 기존 코드와의 호환성을 위해 wide 포맷으로 저장
  ground_truth <- cbind(gs, gc)
  colnames(ground_truth) <- c("gs", "gc")
  write.csv(ground_truth,
            file.path(save_dir, paste0("GT_SI", round(signal_intensity, 2), ".csv")), 
            row.names = FALSE, 
            quote = TRUE)
}



islsimulation <- function(signal_intensity,
                          n = 1,
                          poipar,
                          senna) {
  cisl <- senna
  ref <- cisl@Coord[["Spatial"]]
  csd <- ref[["distance"]]
  csd[csd >= 0] <- 0
  csd <- -csd
  csd <- (csd - min(csd)) / (max(csd) - min(csd))
  grad_index <- c(signal_intensity * poipar[1:500],
                  - signal_intensity * poipar[501:1000],
                  rep(0, 1e3))
  params <- outer(csd, grad_index) + 
    matrix(rep(poipar, nrow(ref)), byrow = TRUE, nrow = nrow(ref)) 

  params[params < 0] <- 0

  save_dir <- file.path(getwd(), "data", "sim.generated", "sim.isl")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  set.seed(seed = n)
  isl_count <- t(
    matrix(
      rpois(length(params), lambda = params), 
      nrow = nrow(params),
      dimnames = list(rownames(ref),
                      colnames(cisl@Gene[["Spatial"]]))))
  
  write.csv(isl_count,
            file.path(save_dir, paste0("sim_isl_SI", round(signal_intensity, 2), ".csv")),
            row.names = TRUE,
            quote = TRUE)

  isl_perm <- t(
    matrix(
      apply(isl_count, 1, sample),
      ncol = nrow(isl_count),
      dimnames = list(colnames(isl_count),
                      rownames(isl_count))))
  
  write.csv(isl_perm,
            file.path(save_dir, paste0("perm_isl_SI", round(signal_intensity, 2), ".csv")),
            row.names = TRUE,
            quote = TRUE)
  
  # Ground truth
  gs <- rownames(isl_count)[1:1000]
  gc <- rownames(isl_count)[1001:2000]
  # 기존 코드와의 호환성을 위해 wide 포맷으로 저장
  ground_truth <- cbind(gs, gc)
  colnames(ground_truth) <- c("gs", "gc")
  write.csv(ground_truth,
            file.path(save_dir, paste0("GT_SI",round(signal_intensity, 2), ".csv")), 
            row.names = FALSE, 
            quote = TRUE)  
}

#buffer

# ==============================
# Raw counts generation
# ==============================

cat("Preprocessing... \n")

t0 <- Sys.time()
s1 <- Load10X_Spatial(
  data.dir = dpath,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Visium",
  filter.matrix = TRUE,
  image = NULL)

if(min(s1$nCount_Spatial) == 0) s1 <- subset(s1, nCount_Spatial > 0)
surt <- s1
s1 <- SCTransform(s1, assay = "Spatial", variable.features.n = 2000, verbose = TRUE)

sen <- SENNA_Visium(s1,
                    slice_name = "Visium")
rm(s1); gc()

#write.csv(sen@Coord[["Spatial"]],
#          paste0(getwd(), "/data/sim.generated/coord.csv"),
#          row.names = TRUE,
#          quote = TRUE)

#AppDat(sen, 
#       image_path = paste0(path, "dataset/SNUH/Thymus9/outs/spatial/"),
#       image_resolution = "lowres")
#knot_picker()

cat("Loading knots... \n")

prsim <- read.csv(paste0(file.path(getwd(), "data", "knots", "svgp.csv")),
                  header = TRUE)
rgsim <- read.csv(paste0(file.path(getwd(), "data", "knots", "svgr.csv")),
                  header = TRUE)
issim <- read.csv(paste0(file.path(getwd(), "data", "knots", "svgi.csv")),
                  header = TRUE)

cat("Preparing SENNA objects... \n")

seni <- TrimmedCurve(sen, issim, type = "islet")
senr <- FullCurve(sen, rgsim, type = "spline")
senp <- FullCurve(sen, prsim, type = "spline")

seni <- GetCurveParam(seni, cores = 1L)
senr <- GetCurveParam(senr, cores = 1L)
senp <- GetCurveParam(senp, cores = 1L)

seni <- TissueRegionation(seni, cores = 1L)
senr <- TissueRegionation(senr, cores = 1L)


cat("Starting simulation... \n")


set.seed(42)
baseline <- rnorm(n = 2000, mean = 10, sd = 1)

n_cores <- parallel::detectCores() - 2
if (is.na(n_cores) || n_cores < 1) n_cores <- 1

cat("Preparing parallel processing...\n" )

{
  cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")
  registerDoParallel(cl)

  clusterExport(cl, varlist = c("progsimulation", "baseline", "senp"), envir = environment())
  clusterEvalQ(cl, { suppressPackageStartupMessages(suppressWarnings(library(SENNA))) })
  }

cat("Generating simulated data for progression pattern... \n")

prog <- foreach(si = seq(from = 0.1, to = 1, length.out = 30)) %dopar% {
  progsimulation(signal_intensity = si,
                 n = 42,
                 poipar = baseline,
                 senna = senp)
}

cat("Prog Finished. Cleaning memory... \n")
stopCluster(cl)
gc()
cat("Memory cleaned. \n")

cat("Preparing parallel processing...\n" )

{
  cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")
  registerDoParallel(cl)

  clusterExport(cl, varlist = c("regiosimulation", "baseline", "senr"), envir = environment())
  clusterEvalQ(cl, { suppressPackageStartupMessages(suppressWarnings(library(SENNA))) })
  }

cat("Generating simulated data for regionation pattern... \n")

regio <- foreach(si = seq(from = 0.1, to = 1, length.out = 30)) %dopar% {
  regiosimulation(signal_intensity = si,
                 n = 42,
                 poipar = baseline,
                 senna = senr)
}

cat("Regio Finished. Cleaning memory... \n")
stopCluster(cl)
gc()
cat("Memory cleaned. \n")


cat("Preparing parallel processing...\n" )

{
  cl <- makeCluster(n_cores, type = "PSOCK", outfile = "")
  registerDoParallel(cl)

  clusterExport(cl, varlist = c("islsimulation", "baseline", "seni"), envir = environment())
  clusterEvalQ(cl, { suppressPackageStartupMessages(suppressWarnings(library(SENNA))) })
  }

cat("Generating simulated data for islet pattern... \n")

islet <- foreach(si = seq(from = 0.1, to = 1, length.out = 30)) %dopar% {
  islsimulation(signal_intensity = si,
                 n = 42,
                 poipar = baseline,
                 senna = seni)
}

cat("Islet Finished. Cleaning memory... \n")
stopCluster(cl)
gc()
cat("Memory cleaned. \n")


cat("Simulation completed. \n")
cat(paste0("Elapsed time: ", Sys.time() - t0, "\n"))