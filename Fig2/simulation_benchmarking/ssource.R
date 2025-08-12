## souce for SENNA simulation

rescalecp <- function(senna){
  coord_mat <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
  fx <- senna@CurveAxis[["fun"]][["x.coef"]]
  fy <- senna@CurveAxis[["fun"]][["y.coef"]]
  knots <- senna@CurveAxis[["knots"]]
  kn <- senna@CurveAxis[["fun"]][["t"]]
  
  if("trimmed" %in% senna@CurveAxis[["type"]]) {
    ll <- lapply(kn[-length(kn)],
                 function(k){
                   qf <- function(x){
                     sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                            (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                   }
                   return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                 })
    ll <- c(0, unlist(ll))
    cll <- cumsum(ll)
    tv <- coord_mat[["t"]]
    tprime <- 
      ifelse(is.nan(tv),
             tv,
             {
               fidx <- floor(tv)
               fidx[fidx == 0] <- 1L
               fidx[fidx > length(cll)] <- length(cll)
               
               (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                 cll[length(cll)]
             })
    senna@Coord[["Spatial"]][["tprime"]] <- tprime
  } else {
    kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
    tv <- coord_mat[["t"]]
    tmin <- min(tv); tmax <- max(tv)
    llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                    ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
    lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                    ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
    ll <- lapply(2:(mk-1),
                 function(k){
                   qf <- function(x){
                     sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                            (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                   }
                   return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                 })
    ll <- c(0, llw, unlist(ll), lhi)
    cll <- cumsum(ll)
    tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
    tv <- ifelse(tv > (max(kn) - 1L), 
                 (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                 tv)
    tv <- tv + 1L
    tprime <- {
      fidx <- floor(tv)
      fidx[fidx == 0] <- 1L
      fidx[fidx >= length(cll)] <- length(cll) - 1L
      (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
        cll[length(cll)]
    }
    senna@Coord[["Spatial"]][["tprime"]] <- tprime
  }
  
  return(senna)
}










progsimulation <- function(signal_intensity,
                           n = 1,
                           poipar,
                           knot,
                           seurat,
                           senna,
                           giotto,
                           sparkarg,
                           pval) {
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
  
  set.seed(seed = n)
  prog_count <- t(
    matrix(
      rpois(length(params), lambda = params), 
      nrow = nrow(params),
      dimnames = list(rownames(ref),
                      colnames(cprog@Gene[["Spatial"]]))))
  prog_perm <- t(
    matrix(
      apply(prog_count, 1, sample),
      ncol = nrow(prog_count),
      dimnames = list(colnames(prog_count),
                      rownames(prog_count))))
  
  # Ground truth
  gs <- rownames(prog_count)[1:1000]
  gc <- rownames(prog_count)[1001:2000]
  
  ## SENNA
  surt_prog <- CreateSeuratObject(
    counts = as.sparse(prog_count),
    assay = "Spatial",
    meta.data = seurat@meta.data)
  surt_prog@images <- seurat@images
  surt_prog@meta.data$nCount_Spatial <- colSums(prog_count)
  surt_prog@meta.data$nFeature_Spatial <- colSums(prog_count > 0)
  surt_prog <- NormalizeData(surt_prog, verbose = FALSE)
  surt_prog <- ScaleData(surt_prog, 
                         features = c(gs, gc),
                         verbose = FALSE)
  sensim <- SENNA_Visium(surt_prog,
                         assay = "Spatial",
                         all_genes = TRUE,
                         slice_name = "Visium")
  sensim <- FullCurve(senna = sensim,
                      knot_df = knot,
                      type = "spline")
  sensim <- GetCovariate(sensim)
  sensim <- ProgSVGs(sensim,
                     weight = "gaussian",
                     FDR_level = pval,
                     grad_cutoff = 0.1)
  sensim <- c(sensim@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]],
              sensim@Gene[["P.SVGs"]][["Variable_gene"]][["negative"]])
  
  surt_perm <- CreateSeuratObject(
    counts = as.sparse(prog_perm),
    assay = "Spatial",
    meta.data = seurat@meta.data)
  surt_perm@images <- seurat@images
  surt_perm@meta.data$nCount_Spatial <- colSums(prog_perm)
  surt_perm@meta.data$nFeature_Spatial <- colSums(prog_perm > 0)
  surt_perm <- NormalizeData(surt_perm, verbose = FALSE)
  surt_perm <- ScaleData(surt_perm, 
                         features = c(gs, gc),
                         verbose = FALSE)
  
  senprm <- SENNA_Visium(surt_perm,
                         assay = "Spatial",
                         all_genes = TRUE,
                         slice_name = "Visium")
  senprm <- FullCurve(senna = senprm,
                      knot_df = knot,
                      type = "spline")
  senprm <- GetCovariate(senprm)
  senprm <- ProgSVGs(senprm,
                     weight = "gaussian",
                     FDR_level = pval,
                     grad_cutoff = 0.1)
  senprm <- c(senprm@Gene[["P.SVGs"]][["Variable_gene"]][["positive"]],
              senprm@Gene[["P.SVGs"]][["Variable_gene"]][["negative"]])
  
  rm(senna); rm(seurat); rm(knot); gc()
  
  
  ## Seurat Moran's I
  morsim <- Seurat::RunMoransI(
    data = surt_prog[["Spatial"]]$scale.data,
    pos = GetTissueCoordinates(object = surt_prog[["Visium"]])[,1:2],
    verbose = FALSE)
  morsim <- rownames(dplyr::filter(morsim, p.value <= pval))
  
  rm(surt_prog); gc()
  
  morprm <- Seurat::RunMoransI(
    data = surt_perm[["Spatial"]]$scale.data,
    pos = GetTissueCoordinates(object = surt_perm[["Visium"]])[,1:2],
    verbose = FALSE)
  morprm <- rownames(dplyr::filter(morprm, p.value <= pval))
  
  rm(surt_perm); gc()
  
  
  ## SPARK-X
  spksim <- sparkx(prog_count,
                   sparkarg,
                   option ="mixture", 
                   verbose = FALSE)
  spksim <- rownames(dplyr::filter(spksim$res_mtest, 
                                   adjustedPval <= pval))
  
  spkprm <- sparkx(prog_perm,
                   sparkarg,
                   option ="mixture", 
                   verbose = FALSE)
  spkprm <- rownames(dplyr::filter(spkprm$res_mtest, 
                                   adjustedPval <= pval))
  
  gc()
  
  
  ## Giotto kmeans
  giksim <- giotto
  giksim@raw_exprs <- as.sparse(prog_count)
  giksim <- normalizeGiotto(gobject = giksim)
  giksim <- createSpatialNetwork(gobject = giksim, minimum_k = 0)
  giksim <- binSpect(giksim, bin_method = 'kmeans', verbose = FALSE)
  giksim <- filter(giksim, adj.p.value <= pval)[["genes"]]
  
  gikprm <- giotto
  gikprm@raw_exprs <- as.sparse(prog_perm)
  gikprm <- normalizeGiotto(gobject = gikprm)
  gikprm <- createSpatialNetwork(gobject = gikprm, minimum_k = 0)
  gikprm <- binSpect(gikprm, bin_method = 'kmeans', verbose = FALSE)
  gikprm <- filter(gikprm, adj.p.value <= pval)[["genes"]]
  
  
  ## Giotto rank
  girsim <- giotto
  girsim@raw_exprs <- as.sparse(prog_count)
  girsim <- normalizeGiotto(gobject = girsim)
  girsim <- createSpatialNetwork(gobject = girsim, minimum_k = 0)
  girsim <- binSpect(girsim, bin_method = 'rank', verbose = FALSE)
  girsim <- filter(girsim, adj.p.value <= pval)[["genes"]]
  
  girprm <- giotto
  girprm@raw_exprs <- as.sparse(prog_perm)
  girprm <- normalizeGiotto(gobject = girprm)
  girprm <- createSpatialNetwork(gobject = girprm, minimum_k = 0)
  girprm <- binSpect(girprm, bin_method = 'rank', verbose = FALSE)
  girprm <- filter(girprm, adj.p.value <= pval)[["genes"]]
  
  rm(giotto); gc()
  
  ## Output
  sim_res <- tibble(
    ID = c("DR_sim", "FPR_sim", "DR_prm"),
    SENNA = c(length(intersect(gs, sensim)) / length(gs),
              length(intersect(gc, sensim)) / length(gc),
              length(intersect(gs, senprm)) / length(gs)),
    Seurat = c(length(intersect(gs, morsim)) / length(gs),
               length(intersect(gc, morsim)) / length(gc),
               length(intersect(gs, morprm)) / length(gs)),
    SPARK_X = c(length(intersect(gs, spksim)) / length(gs),
                length(intersect(gc, spksim)) / length(gc),
                length(intersect(gs, spkprm)) / length(gs)),
    Giotto_kmeans = c(length(intersect(gs, giksim)) / length(gs),
                      length(intersect(gc, giksim)) / length(gc),
                      length(intersect(gs, gikprm)) / length(gs)),
    Giotto_rank = c(length(intersect(gs, girsim)) / length(gs),
                    length(intersect(gc, girsim)) / length(gc),
                    length(intersect(gs, girprm)) / length(gs)),
    SI = rep(signal_intensity, 3)
  )
  
  return(sim_res)
}










regiosimulation <- function(signal_intensity,
                            n = 1,
                            poipar,
                            knot,
                            seurat,
                            senna,
                            giotto,
                            sparkarg,
                            pval) {
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
  
  set.seed(seed = n)
  regio_count <- t(
    matrix(
      rpois(length(params), lambda = params), 
      nrow = nrow(params),
      dimnames = list(rownames(ref),
                      colnames(cregio@Gene[["Spatial"]]))))
  
  regio_perm <- t(
    matrix(
      apply(regio_count, 1, sample),
      ncol = nrow(regio_count),
      dimnames = list(colnames(regio_count),
                      rownames(regio_count))))
  
  
  # Ground truth
  gs <- rownames(regio_count)[1:1000]
  gc <- rownames(regio_count)[1001:2000]
  
  ## SENNA
  surt_regio <- CreateSeuratObject(
    counts = as.sparse(regio_count),
    assay = "Spatial",
    meta.data = seurat@meta.data)
  surt_regio@images <- seurat@images
  surt_regio@meta.data$nCount_Spatial <- colSums(regio_count)
  surt_regio@meta.data$nFeature_Spatial <- colSums(regio_count > 0)
  surt_regio <- NormalizeData(surt_regio, verbose = FALSE)
  surt_regio <- ScaleData(surt_regio, 
                          features = c(gs, gc),
                          verbose = FALSE)
  sensim <- SENNA_Visium(surt_regio,
                         assay = "Spatial",
                         all_genes = TRUE,
                         slice_name = "Visium")
  sensim <- FullCurve(senna = sensim,
                      knot_df = knot,
                      type = "spline")
  sensim <- GetCovariate(sensim)
  sensim <- TissueRegionation(sensim)
  sensim <- RegionSVGs(sensim,
                       FDR_level = pval,
                       grad_cutoff = 0.1)
  sensim <- c(sensim@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
              sensim@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])
  
  surt_perm <- CreateSeuratObject(
    counts = as.sparse(regio_perm),
    assay = "Spatial",
    meta.data = seurat@meta.data)
  surt_perm@images <- seurat@images
  surt_perm@meta.data$nCount_Spatial <- colSums(regio_perm)
  surt_perm@meta.data$nFeature_Spatial <- colSums(regio_perm > 0)
  DefaultAssay(surt_perm) <- "Spatial"
  surt_perm <- NormalizeData(surt_perm, verbose = FALSE)
  surt_perm <- ScaleData(surt_perm, 
                         features = c(gs, gc),
                         verbose = FALSE)
  surt_perm <- FindVariableFeatures(surt_perm, 
                                    variable.features.n = 2000,
                                    verbose = FALSE)
  senprm <- SENNA_Visium(surt_perm,
                         assay = "Spatial",
                         all_genes = TRUE,
                         slice_name = "Visium")
  senprm <- FullCurve(senna = senprm,
                      knot_df = knot,
                      type = "spline")
  senprm <- GetCovariate(senprm)
  senprm <- TissueRegionation(senprm)
  senprm <- RegionSVGs(senprm,
                       FDR_level = pval,
                       grad_cutoff = 0.1)
  senprm <- c(senprm@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
              senprm@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])
  
  rm(senna); rm(seurat); rm(knot); gc()
  
  ## Seurat Moran's I
  morsim <- Seurat::RunMoransI(
    data = surt_regio[["Spatial"]]$scale.data,
    pos = GetTissueCoordinates(object = surt_regio[["Visium"]])[,1:2],
    verbose = FALSE)
  morsim <- rownames(dplyr::filter(morsim, p.value <= pval))
  
  rm(surt_regio); gc()
  
  morprm <- Seurat::RunMoransI(
    data = surt_perm[["Spatial"]]$scale.data,
    pos = GetTissueCoordinates(object = surt_perm[["Visium"]])[,1:2],
    verbose = FALSE)
  morprm <- rownames(dplyr::filter(morprm, p.value <= pval))
  
  rm(surt_perm); gc()
  
  
  ## SPARK-X
  spksim <- sparkx(regio_count,
                   sparkarg,
                   option ="mixture", 
                   verbose = FALSE)
  spksim <- rownames(dplyr::filter(spksim$res_mtest, 
                                   adjustedPval <= pval))
  
  spkprm <- sparkx(regio_perm,
                   sparkarg,
                   option ="mixture", 
                   verbose = FALSE)
  spkprm <- rownames(dplyr::filter(spkprm$res_mtest, 
                                   adjustedPval <= pval))
  
  gc()
  
  
  ## Giotto kmeans
  giksim <- giotto
  giksim@raw_exprs <- as.sparse(regio_count)
  giksim <- normalizeGiotto(gobject = giksim)
  giksim <- createSpatialNetwork(gobject = giksim, minimum_k = 0)
  giksim <- binSpect(giksim, bin_method = 'kmeans', verbose = FALSE)
  giksim <- filter(giksim, adj.p.value <= pval)[["genes"]]
  
  gikprm <- giotto
  gikprm@raw_exprs <- as.sparse(regio_perm)
  gikprm <- normalizeGiotto(gobject = gikprm)
  gikprm <- createSpatialNetwork(gobject = gikprm, minimum_k = 0)
  gikprm <- binSpect(gikprm, bin_method = 'kmeans', verbose = FALSE)
  gikprm <- filter(gikprm, adj.p.value <= pval)[["genes"]]
  
  
  ## Giotto rank
  girsim <- giotto
  girsim@raw_exprs <- as.sparse(regio_count)
  girsim <- normalizeGiotto(gobject = girsim)
  girsim <- createSpatialNetwork(gobject = girsim, minimum_k = 0)
  girsim <- binSpect(girsim, bin_method = 'rank', verbose = FALSE)
  girsim <- filter(girsim, adj.p.value <= pval)[["genes"]]
  
  girprm <- giotto
  girprm@raw_exprs <- as.sparse(regio_perm)
  girprm <- normalizeGiotto(gobject = girprm)
  girprm <- createSpatialNetwork(gobject = girprm, minimum_k = 0)
  girprm <- binSpect(girprm, bin_method = 'rank', verbose = FALSE)
  girprm <- filter(girprm, adj.p.value <= pval)[["genes"]]
  
  rm(giotto); gc()
  
  ## Output
  sim_res <- tibble(
    ID = c("DR_sim", "FPR_sim", "DR_prm"),
    SENNA = c(length(intersect(gs, sensim)) / length(gs),
              length(intersect(gc, sensim)) / length(gc),
              length(intersect(gs, senprm)) / length(gs)),
    Seurat = c(length(intersect(gs, morsim)) / length(gs),
               length(intersect(gc, morsim)) / length(gc),
               length(intersect(gs, morprm)) / length(gs)),
    SPARK_X = c(length(intersect(gs, spksim)) / length(gs),
                length(intersect(gc, spksim)) / length(gc),
                length(intersect(gs, spkprm)) / length(gs)),
    Giotto_kmeans = c(length(intersect(gs, giksim)) / length(gs),
                      length(intersect(gc, giksim)) / length(gc),
                      length(intersect(gs, gikprm)) / length(gs)),
    Giotto_rank = c(length(intersect(gs, girsim)) / length(gs),
                    length(intersect(gc, girsim)) / length(gc),
                    length(intersect(gs, girprm)) / length(gs)),
    SI = rep(signal_intensity, 3)
  )
  
  return(sim_res)
}










islsimulation <- function(signal_intensity,
                          n = 1,
                          poipar,
                          knot,
                          seurat,
                          senna,
                          giotto,
                          sparkarg,
                          pval) {
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
  
  set.seed(seed = n)
  isl_count <- t(
    matrix(
      rpois(length(params), lambda = params), 
      nrow = nrow(params),
      dimnames = list(rownames(ref),
                      colnames(cisl@Gene[["Spatial"]]))))
  isl_perm <- t(
    matrix(
      apply(isl_count, 1, sample),
      ncol = nrow(isl_count),
      dimnames = list(colnames(isl_count),
                      rownames(isl_count))))
  
  
  # Ground truth
  gs <- rownames(isl_count)[1:1000]
  gc <- rownames(isl_count)[1001:2000]
  
  ## SENNA
  surt_isl <- CreateSeuratObject(
    counts = as.sparse(isl_count),
    assay = "Spatial",
    meta.data = seurat@meta.data)
  surt_isl@images <- seurat@images
  surt_isl@meta.data$nCount_Spatial <- colSums(isl_count)
  surt_isl@meta.data$nFeature_Spatial <- colSums(isl_count > 0)
  surt_isl <- NormalizeData(surt_isl, verbose = FALSE)
  surt_isl <- ScaleData(surt_isl, 
                        features = c(gs, gc),
                        verbose = FALSE)
  sensim <- SENNA_Visium(surt_isl,
                         assay = "Spatial",
                         all_genes = TRUE,
                         slice_name = "Visium")
  sensim <- TrimmedCurve(senna = sensim,
                         knot_df = knot,
                         type = "islet")
  sensim <- GetCovariate(sensim)
  sensim <- TissueRegionation(sensim)
  sensim <- RegionSVGs(sensim,
                       direction = -1,
                       FDR_level = pval,
                       grad_cutoff = 0.1)
  sensim <- c(sensim@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
              sensim@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])
  
  surt_perm <- CreateSeuratObject(
    counts = as.sparse(isl_perm),
    assay = "Spatial",
    meta.data = seurat@meta.data)
  surt_perm@images <- seurat@images
  surt_perm@meta.data$nCount_Spatial <- colSums(isl_perm)
  surt_perm@meta.data$nFeature_Spatial <- colSums(isl_perm > 0)
  surt_perm <- NormalizeData(surt_perm, verbose = FALSE)
  surt_perm <- ScaleData(surt_perm, 
                         features = c(gs, gc),
                         verbose = FALSE)
  senprm <- SENNA_Visium(surt_perm,
                         assay = "Spatial",
                         all_genes = TRUE,
                         slice_name = "Visium")
  senprm <- TrimmedCurve(senna = senprm,
                         knot_df = knot,
                         type = "islet")
  senprm <- GetCovariate(senprm)
  senprm <- TissueRegionation(senprm)
  senprm <- RegionSVGs(senprm,
                       direction = -1,
                       FDR_level = pval,
                       grad_cutoff = 0.1)
  senprm <- c(senprm@Gene[["R.SVGs"]][["Variable_gene"]][["positive"]],
              senprm@Gene[["R.SVGs"]][["Variable_gene"]][["negative"]])
  
  rm(senna); rm(seurat); rm(knot); gc()
  
  ## Seurat Moran's I
  morsim <- Seurat::RunMoransI(
    data = surt_isl[["Spatial"]]$scale.data,
    pos = GetTissueCoordinates(object = surt_isl[["Visium"]])[,1:2],
    verbose = FALSE)
  morsim <- rownames(dplyr::filter(morsim, p.value <= pval))
  
  rm(srut_isl); gc()
  
  morprm <- Seurat::RunMoransI(
    data = surt_perm[["Spatial"]]$scale.data,
    pos = GetTissueCoordinates(object = surt_perm[["Visium"]])[,1:2],
    verbose = FALSE)
  morprm <- rownames(dplyr::filter(morprm, p.value <= pval))
  
  rm(surt_perm); gc()
  
  
  ## SPARK-X
  spksim <- sparkx(isl_count,
                   sparkarg,
                   option ="mixture", 
                   verbose = FALSE)
  spksim <- rownames(dplyr::filter(spksim$res_mtest, 
                                   adjustedPval <= pval))
  
  spkprm <- sparkx(isl_perm,
                   sparkarg,
                   option ="mixture", 
                   verbose = FALSE)
  spkprm <- rownames(dplyr::filter(spkprm$res_mtest, 
                                   adjustedPval <= pval))
  
  gc()
  
  
  ## Giotto kmeans
  giksim <- giotto
  giksim@raw_exprs <- as.sparse(isl_count)
  giksim <- normalizeGiotto(gobject = giksim)
  giksim <- createSpatialNetwork(gobject = giksim, minimum_k = 0)
  giksim <- binSpect(giksim, bin_method = 'kmeans', verbose = FALSE)
  giksim <- filter(giksim, adj.p.value <= pval)[["genes"]]
  
  gikprm <- giotto
  gikprm@raw_exprs <- as.sparse(isl_perm)
  gikprm <- normalizeGiotto(gobject = gikprm)
  gikprm <- createSpatialNetwork(gobject = gikprm, minimum_k = 0)
  gikprm <- binSpect(gikprm, bin_method = 'kmeans', verbose = FALSE)
  gikprm <- filter(gikprm, adj.p.value <= pval)[["genes"]]
  
  
  ## Giotto rank
  girsim <- giotto
  girsim@raw_exprs <- as.sparse(isl_count)
  girsim <- normalizeGiotto(gobject = girsim)
  girsim <- createSpatialNetwork(gobject = girsim, minimum_k = 0)
  girsim <- binSpect(girsim, bin_method = 'rank', verbose = FALSE)
  girsim <- filter(girsim, adj.p.value <= pval)[["genes"]]
  
  girprm <- giotto
  girprm@raw_exprs <- as.sparse(isl_perm)
  girprm <- normalizeGiotto(gobject = girprm)
  girprm <- createSpatialNetwork(gobject = girprm, minimum_k = 0)
  girprm <- binSpect(girprm, bin_method = 'rank', verbose = FALSE)
  girprm <- filter(girprm, adj.p.value <= pval)[["genes"]]
  
  rm(giotto); gc()
  
  ## Output
  sim_res <- tibble(
    ID = c("DR_sim", "FPR_sim", "DR_prm"),
    SENNA = c(length(intersect(gs, sensim)) / length(gs),
              length(intersect(gc, sensim)) / length(gc),
              length(intersect(gs, senprm)) / length(gs)),
    Seurat = c(length(intersect(gs, morsim)) / length(gs),
               length(intersect(gc, morsim)) / length(gc),
               length(intersect(gs, morprm)) / length(gs)),
    SPARK_X = c(length(intersect(gs, spksim)) / length(gs),
                length(intersect(gc, spksim)) / length(gc),
                length(intersect(gs, spkprm)) / length(gs)),
    Giotto_kmeans = c(length(intersect(gs, giksim)) / length(gs),
                      length(intersect(gc, giksim)) / length(gc),
                      length(intersect(gs, gikprm)) / length(gs)),
    Giotto_rank = c(length(intersect(gs, girsim)) / length(gs),
                    length(intersect(gc, girsim)) / length(gc),
                    length(intersect(gs, girprm)) / length(gs)),
    SI = rep(signal_intensity, 3)
  )
  
  return(sim_res)
}





