
# This script contains amended functions from the packages APE (parafit), PACo and Mantel test to detect cophylogenetic and phylogenetic signals, by taking into account biogeography or not


#### Mantel test ####

mantel_test <- function(network, tree_A, tree_B = NULL, method="Jaccard_binary", 
                        nperm=1000, correlation="Pearson", permutation="shuffle",
                        verbose=TRUE, table_biogeo=NULL){
  
  if (!correlation %in% c("Pearson", "Spearman")) {stop("\"correlation\" must be among 'Pearson' or 'Spearman'.")}
  if (!is.numeric(nperm)) {stop("Please provide a numeric number of permutations (\"nperm\").")}
  
  if (!inherits(tree_A, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!is.null(tree_B)) {if (!inherits(tree_B, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}}
  
  if (is.null(method)) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', and 'UniFrac_unweighted'.")}
  if (method %in% c("GUniFrac", "UniFrac_unweighted", "PBLM", "PBLM_binary")) {if (is.null(tree_B)) stop("Please provide a phylogenetic tree \"tree_B\" for guild B.")}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "Bray-Curtis", "GUniFrac", "UniFrac_unweighted")) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', and 'UniFrac_unweighted'.")}
  
  if (nrow(network)<2){stop("Please provide a \"network\" with at least 2 species in clade B.")}
  if (ncol(network)<2){stop("Please provide a \"network\" with at least 2 species in clade A.")}
  
  if (!permutation %in% c("shuffle","nbpartners", "biogeo")) {stop("Please provide a type of \"permutation\" among 'shuffle' and 'nbpartners'.")}
  
  if (permutation == "biogeo") { 
    if (!all(table_biogeo$species %in% colnames(network))){
      stop("Provide a table_biogeo with all the species from tree_A.")
    }
    table_biogeo <- table_biogeo[match(colnames(network), table_biogeo$species), ]
  }
  
  
  # Only keep species with at least 1 interaction
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  
  # A in columns and B in rows
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  names(nb_A) <- "nb_A"
  names(nb_B) <- "nb_B"
  
  # Check names
  if (all(is.null(colnames(network)))|all(is.null(rownames(network)))) {stop("Please provide a \"network\" with row names and columns names matching the species names.")}
  
  if (!all(colnames(network) %in% tree_A$tip.label)){stop("Please provide a \"tree_A\" for all the species in clade A (the columns of the intercation network).")}
  
  tree_A <- ape::drop.tip(tree_A,tip=tree_A$tip.label[which(!tree_A$tip.label %in% colnames(network))])
  
  if (!is.rooted(tree_A)){tree_A <- phytools::midpoint.root(tree_A) }
  
  network <- network[1:nrow(network),tree_A$tip.label]
  
  compute_eco_dist <- function(network, method){
    # binary Jaccard distances
    if (method=="Jaccard_binary"){
      jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=TRUE))
      eco_A <- jaccard_A
    }
    
    # quantitative Jaccard distances
    if (method=="Jaccard_weighted"){
      jaccard_A <- as.matrix(vegan::vegdist(t(network), "jaccard", binary=FALSE))
      eco_A <- jaccard_A
    }
    
    # Bray-Curtis dissimilarity 
    if (method=="Bray-Curtis"){
      bray_A <- as.matrix(vegan::vegdist(t(network), "bray", binary=FALSE))
      eco_A <- bray_A
    }
    
    # Unifrac (generalized UniFrac, with alpha=0.5)
    if (method=="GUniFrac"){
      unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
      index=1
      eco_A <- unifrac_A$unifracs[,,index]
    }
    
    # Unifrac (unweighted UniFrac)
    if (method=="UniFrac_unweighted"){
      unifrac_A <- GUniFrac::GUniFrac(t(network), tree = tree_B, alpha=c(0.5))
      index=2
      eco_A <- unifrac_A$unifracs[,,index]
    }
    
    # Degree
    if (method=="degree"){
      network_binary <- network
      network_binary[network_binary>0] <- 1
      eco_A <- as.matrix(dist(colSums(network_binary)))
    }
    return(eco_A)
  }
  
  eco_A <- compute_eco_dist(network, method)
  
  # Perform Mantel test:
  
  # cophenetic distances
  cophe_A <- ape::cophenetic.phylo(tree_A)
  
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  
  results <- list(R = NA, p_value_upper = NA, p_value_lower = NA, random_correlations = NA)
  
  if (length(unique(as.vector(cophe_A)))<3) {
    if (verbose) print("The phylogenetic distance matrix is composed of only 2 different values (because of polytomies?).")
    return(results)}
  if (length(unique(as.vector(eco_A)))<3) {
    if (verbose) print("The ecological distance matrix is composed of only 1 value (identical patterns of interactions across species?).")
    return(results)}
  
  
  if (correlation=="Pearson") {correlation="pearson"}
  if (correlation=="Spearman") {correlation="spearman"}
  
  original_correlation <- cor(as.vector(as.dist(eco_A)), as.vector(as.dist(cophe_A)), use = "everything", method = correlation)
  
  
  # Make randomizations
  
  random_correlation <- vector(mode = "numeric", length = nperm)
  vector_cophe_A <- as.vector(as.dist(cophe_A))
  
  if (permutation=="shuffle"){
    for (i in 1:nperm){
      random <- sample(1:nb_A)
      eco_A_rand <- eco_A[random,random]
      random_correlation[i] <- cor(as.vector(as.dist(eco_A_rand)), vector_cophe_A, use = "everything", method = correlation)
    }
  }
  
  if (permutation=="nbpartners"){
    for (i in 1:nperm){
      rand_network <- network
      for (k in 1:nb_A){
        rand_network[,k] <- sample(rand_network[,k])
      }
      eco_A_rand <- compute_eco_dist(rand_network)
      random_correlation[i] <- cor(as.vector(as.dist(eco_A_rand)), vector_cophe_A, use = "everything", method = correlation)
    }
  }
  
  if (permutation=="biogeo"){
    for (i in 1:nperm){
      rand_network <- network
      for (zone in unique(table_biogeo$biogeo)){
        w <- which(table_biogeo$biogeo==zone)
        s <- as.numeric(sample(as.character(w))) 
        rand_network[,w] <- as.matrix(rand_network[,s])
      }
      eco_A_rand <- compute_eco_dist(rand_network)
      random_correlation[i] <- cor(as.vector(as.dist(eco_A_rand)), vector_cophe_A, use = "everything", method = correlation)
    }
  }
  
  out <- list(R = original_correlation, 
              p_value_upper = min(c(length(which(c(random_correlation,original_correlation)>=original_correlation))/nperm,1)), 
              p_value_lower = min(c(length(which(c(random_correlation,original_correlation)<=original_correlation))/nperm, 1)),
              random_correlations=random_correlation)
  
  return(out)
}






#### Phylosignal network ####

phylosignal_network <- function(network, tree_A, tree_B=NULL, method = "Jaccard_weighted", nperm = 10000, correlation = "Pearson", only_A = FALSE, 
                                permutation ="shuffle", table_biogeo=NULL){
  
  
  if (is.null(tree_B)) {only_A <- TRUE} 
  
  if (!inherits(tree_A, "phylo")) {stop("object \"tree_A\" is not of class \"phylo\".")}
  if (!is.null(tree_B)) {if (!inherits(tree_B, "phylo")) {stop("object \"tree_B\" is not of class \"phylo\".")}}
  
  if (is.null(method)) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', 'UniFrac_unweighted', 'PBLM', 'PBLM_binary', and 'degree'.")}
  if (method %in% c("GUniFrac", "UniFrac_unweighted")) {if (is.null(tree_B)) stop("Please provide a phylogenetic tree \"tree_B\" for guild B.")}
  if (!method %in% c("Jaccard_weighted","Jaccard_binary", "Bray-Curtis", "GUniFrac", "UniFrac_unweighted", "degree")) {stop("Please provide a \"method\" to compute phylogenetic signals among 'Jaccard_weighted', 'Jaccard_binary', 'Bray-Curtis', 'GUniFrac', 'UniFrac_unweighted', and 'degree'.")}
  
  
  if (!permutation %in% c("shuffle","nbpartners", "biogeo")) {stop("Please provide a type of \"permutation\" among 'shuffle' and 'nbpartners'.")}
  if (permutation!="shuffle") {if (method %in% c("degree")) stop("The argument \"permutation\" is not used for this method.")}
  
  
  if (permutation == "biogeo") { 
    if (!all(table_biogeo$species %in% colnames(network))){
      stop("Provide a table_biogeo with all the species from tree_A.")
    }
    table_biogeo <- table_biogeo[match(colnames(network), table_biogeo$species), ]
  }
  
  
  if (!correlation %in% c("Pearson", "Spearman", "Kendall")) {stop("Please pick a \"correlation\" among Pearson, Spearman, and Kendall.")}
  
  if (nrow(network)<2){stop("Please provide a \"network\" with at least 2 species in clade B.")}
  if (ncol(network)<2){stop("Please provide a \"network\" with at least 2 species in clade A.")}
  
  
  # Only keep species with at least 1 interaction
  network <- network[rowSums(network)>0,]
  network <- network[,colSums(network)>0]
  
  # A in columns and B in rows
  nb_A <- ncol(network)
  nb_B <- nrow(network)
  names(nb_A) <- "nb_A"
  names(nb_B) <- "nb_B"
  
  # Check names
  if (all(is.null(colnames(network)))|all(is.null(rownames(network)))) {stop("Please provide a \"network\" with row names and columns names matching the species names.")}
  
  if (!all(colnames(network) %in% tree_A$tip.label)){stop("Please provide a \"tree_A\" for all the species in clade A (the columns of the interaction network).")}
  if (only_A==FALSE) { if (!all(rownames(network) %in% tree_B$tip.label)){stop("Please provide a \"tree_B\" for all the species in clade B (the rows of the interaction network).")}}
  
  tree_A <- ape::drop.tip(tree_A,tip=tree_A$tip.label[which(!tree_A$tip.label %in% colnames(network))])
  if (only_A==FALSE) { tree_B <- ape::drop.tip(tree_B,tip=tree_B$tip.label[which(!tree_B$tip.label %in% rownames(network))])}
  
  
  if (!is.rooted(tree_A)){tree_A <- phytools::midpoint.root(tree_A) }
  if (only_A==FALSE) { if (!is.rooted(tree_B)){tree_B <- phytools::midpoint.root(tree_B) }}
  
  if (only_A==TRUE) { 
    network <- network[1:nrow(network),tree_A$tip.label]
  } else {
    network <- network[tree_B$tip.label,tree_A$tip.label]
  }
  
  # Mantel tests
  
  if (permutation=="shuffle"){
    mantel_A <- mantel_test(network, tree_A, tree_B, method, nperm, correlation, permutation)
    if (only_A==FALSE) {mantel_B <- mantel_test(t(network), tree_B, tree_A, method, nperm, correlation, permutation)
    }else{mantel_B <- list(R=NA, p_value_upper=NA, p_value_lower=NA)}
  }
  
  
  if (permutation=="nbpartners"){
    mantel_A <- mantel_test(network, tree_A, tree_B, method, nperm, correlation, permutation)
    if (only_A==FALSE) {mantel_B <- mantel_test(t(network), tree_B, tree_A, method, nperm, correlation, permutation)
    }else{mantel_B <- list(R=NA, p_value_upper=NA, p_value_lower=NA)}
  }
  
  if (permutation=="biogeo"){
    mantel_A <- mantel_test(network=network, tree_A=tree_A, tree_B=tree_B, method=method, nperm = nperm, correlation= correlation, table_biogeo=table_biogeo)
    mantel_B <- list(R=NA, p_value_upper=NA, p_value_lower=NA)
  }
  
  out <- list(nb_A=as.integer(nb_A), nb_B=as.integer(nb_B), 
              mantel_cor_A=mantel_A$R, p_upper_A=mantel_A$p_value_upper, p_lower_A=mantel_A$p_value_lower,mantel_cor_A_perm=mantel_A$random_correlation, 
              mantel_cor_B=mantel_B$R, p_upper_B=mantel_B$p_value_upper, p_lower_B=mantel_B$p_value_lower,mantel_cor_B_perm=mantel_B$random_correlation)
  
  return(out)
}



parafit_test <- function (dist_A, dist_B, network, nperm = 999,
                          seed = NULL, correction = "none", silent = TRUE,
                          nullmodel=2, table_biogeo=NULL) {
  epsilon <- sqrt(.Machine$double.eps)
  
  if (is.null(seed)) {
    runif(1)
    seed <- .Random.seed[trunc(runif(1, 1, 626))]
  }
  
  HP <- t(as.matrix(network))
  host.D <- dist_A
  para.D <- dist_B
  # now tree_A on rows and tree_B on columns
  
  if (nullmodel == 3) { 
    if (!all(table_biogeo$species %in% rownames(host.D))){
      stop("Provide a table_biogeo with all the species from tree_A.")
    }
    table_biogeo <- table_biogeo[match(rownames(host.D), table_biogeo$species), ]
  }
  
  
  
  host.D <- as.matrix(host.D)
  host.pc <- pcoa_test(host.D, correction = correction)
  if (host.pc$correction[2] == 1) {
    if (min(host.pc$values[, 2]) < -epsilon) 
      stop("Tree A matrix has negative eigenvalues. Rerun with correction=\"lingoes\" or correction=\"cailliez\"")
    sum.host.values.sq <- sum(host.pc$values[, 1]^2)
    host.vectors <- host.pc$vectors
  } else {
    sum.host.values.sq <- sum(host.pc$values[, 2]^2)
    host.vectors <- host.pc$vectors.cor
  }
  if (!nullmodel %in% c(1,2,3)) {
    stop("Pick a null model among \"1\" or \"2\" or \"3\".")
  }
  n.host <- nrow(host.D)
  para.D <- as.matrix(para.D)
  para.pc <- pcoa_test(para.D, correction = correction)
  if (para.pc$correction[2] == 1) {
    if (min(para.pc$values[, 2]) < -epsilon) 
      stop("Tree B matrix has negative eigenvalues. Rerun with correction=\"lingoes\" or correction=\"cailliez\"")
    sum.para.values.sq <- sum(para.pc$values[, 1]^2)
    para.vectors <- para.pc$vectors
  }else {
    sum.para.values.sq <- sum(para.pc$values[, 2]^2)
    para.vectors <- para.pc$vectors.cor
  }
  n.para <- nrow(para.D)
  if (!silent) 
    cat("n.tree_A =", n.host, ", n.tree_B =", n.para, "\n")
  
  tracemax <- max(sum.host.values.sq, sum.para.values.sq)
  if (n.host == n.para) {
    if (!silent) 
      cat("The function cannot check if matrix network has been entered in the right way.", 
          "\n")
    if (!silent) 
      cat("It will assume that the rows of network correspond to the clade B.", 
          "\n")
  }else {
    temp <- dim(HP)
    if (temp[1] == n.host) {
      if (temp[2] != n.para) 
        stop("Matrices dist_A, dist_B and network not comformable")
    }else if (temp[2] == n.host) {
      if (temp[1] != n.para) 
        stop("Matrices dist_A, dist_B and network not comformable")
      !
        if (!silent) 
          cat("Matrix network has been transposed for comformity with dist_A and dist_B.", 
              "\n")
    }else {
      stop("Matrices dist_A, dist_B and network not comformable")
    }
  }
  p.per.h <- apply(HP, 1, sum)
  h.per.p <- apply(HP, 2, sum)
  mat.4 <- t(host.vectors) %*% HP %*% para.vectors
  global <- sum(mat.4^2)
  if (nperm > 0) {
    set.seed(seed)
    nGT <- 1
    global.perm <- NA
    for (i in 1:nperm) {
      
      if (nullmodel==1) {
        HP.perm <- apply(HP, 2, sample)
      }
      if (nullmodel==2) {
        HP.perm <- HP[sample(1:nrow(HP)),]
      }
      if (nullmodel==3){ # Takes into account biogeography
        HP.perm <- HP
        for (zone in unique(table_biogeo$biogeo)){
          w <- which(table_biogeo$biogeo==zone)
          s <- as.numeric(sample(as.character(w))) 
          HP.perm[w,] <- as.matrix(HP.perm[s,])
        }
        rownames(HP.perm) = rownames(HP)
      }
      
      mat.4.perm <- t(host.vectors) %*% HP.perm %*% para.vectors
      global.perm <- c(global.perm, sum(mat.4.perm^2))
      if (global.perm[i + 1] >= global) 
        nGT <- nGT + 1
    }
    global.perm <- global.perm[-1] # under H0 
    p.global <- nGT/(nperm) # nGT/(nperm + 1)
  }else {
    p.global <- NA
  }
  
  out <- list(ParaFitGlobal = global, 
              p = p.global, 
              ParaFitGlobal_perm=global.perm,
              nperm = nperm,
              nullmodel = nullmodel)
  return(out)
}


pcoa_test <-  function (D, correction = "none", rn = NULL) {
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  bstick.def <- function(n, tot.var = 1, ...) {
    res <- rev(cumsum(tot.var/n:1)/n)
    names(res) <- paste("Stick", seq(len = n), sep = "")
    return(res)
  }
  D <- as.matrix(D)
  n <- nrow(D)
  epsilon <- sqrt(.Machine$double.eps)
  if (length(rn) != 0) {
    names <- rn
  } else {
    names <- rownames(D)
  }
  CORRECTIONS <- c("none", "lingoes", "cailliez")
  correct <- pmatch(correction, CORRECTIONS)
  if (is.na(correct)) 
    stop("Invalid correction method")
  delta1 <- centre((-0.5 * D^2), n)
  trace <- sum(diag(delta1))
  D.eig <- eigen(delta1)
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  D.eig$values[zero.eig] <- 0
  
  if (min.eig > -epsilon) {
    correct <- 1
    eig <- D.eig$values
    k <- length(which(eig > epsilon))
    rel.eig <- eig[1:k]/trace
    cum.eig <- cumsum(rel.eig)
    vectors <- sweep(D.eig$vectors[, 1:k, drop=F], 2, sqrt(eig[1:k]), FUN = "*")
    bs <- bstick.def(k)
    cum.bs <- cumsum(bs)
    res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
    colnames(res) <- c("Eigenvalues", "Relative_eig", "Broken_stick", 
                       "Cumul_eig", "Cumul_br_stick")
    rownames(res) <- 1:nrow(res)
    rownames(vectors) <- names
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                  prefix = "Axis.")
    note <- paste("There were no negative eigenvalues. No correction was applied")
    out <- (list(correction = c(correction, correct), note = note, 
                 values = res, vectors = vectors, trace = trace))
  }else {
    k <- n
    eig <- D.eig$values
    rel.eig <- eig/trace
    rel.eig.cor <- (eig - min.eig)/(trace - (n - 1) * min.eig)
    if (length(zero.eig)) 
      rel.eig.cor <- c(rel.eig.cor[-zero.eig[1]], 0)
    cum.eig.cor <- cumsum(rel.eig.cor)
    k2 <- length(which(eig > epsilon))
    k3 <- length(which(rel.eig.cor > epsilon))
    vectors <- sweep(D.eig$vectors[, 1:k2, drop=F], 2, sqrt(eig[1:k2]), FUN = "*")
    if ((correct == 2) | (correct == 3)) {
      if (correct == 2) {
        c1 <- -min.eig
        note <- paste("Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 -", 
                      c1, ", except diagonal elements")
        D <- -0.5 * (D^2 + 2 * c1)
      }else if (correct == 3) {
        delta2 <- centre((-0.5 * D), n)
        upper <- cbind(matrix(0, n, n), 2 * delta1)
        lower <- cbind(-diag(n), -4 * delta2)
        sp.matrix <- rbind(upper, lower)
        c2 <- max(Re(eigen(sp.matrix, symmetric = FALSE, 
                           only.values = TRUE)$values))
        note <- paste("Cailliez correction applied to negative eigenvalues: D' = -0.5*(D +", 
                      c2, ")^2, except diagonal elements")
        D <- -0.5 * (D + c2)^2
      }
      diag(D) <- 0
      mat.cor <- centre(D, n)
      toto.cor <- eigen(mat.cor)
      trace.cor <- sum(diag(mat.cor))
      min.eig.cor <- min(Re(toto.cor$values))
      zero.eig.cor <- which((Re(toto.cor$values) < epsilon) & 
                              (Re(toto.cor$values) > -epsilon)) # idem
      toto.cor$values[zero.eig.cor] <- 0
      toto.cor$values <- Re(toto.cor$values)
      toto.cor$vectors <- Re(toto.cor$vectors) # idem
      
      if (min.eig.cor > -epsilon) {
        eig.cor <- toto.cor$values
        rel.eig.cor <- eig.cor[1:k]/trace.cor
        cum.eig.cor <- cumsum(rel.eig.cor)
        k2 <- length(which(eig.cor > epsilon))
        vectors.cor <- sweep(toto.cor$vectors[, 1:k2,drop=F],2, sqrt(eig.cor[1:k2]), FUN = "*")
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
        bs <- bstick.def(k2)
        bs <- c(bs, rep(0, (k - k2)))
        cum.bs <- cumsum(bs)
      }else {
        if (correct == 2) 
          cat("Problem! Negative eigenvalues are still present after Lingoes", 
              "\n")
        if (correct == 3) 
          cat("Problem! Negative eigenvalues are still present after Cailliez", 
              "\n")
        rel.eig.cor <- cum.eig.cor <- bs <- cum.bs <- rep(NA, 
                                                          n)
        vectors.cor <- matrix(NA, n, 2)
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
      }
      res <- data.frame(eig[1:k], eig.cor[1:k], rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Corr_eig", "Rel_corr_eig", 
                         "Broken_stick", "Cum_corr_eig", "Cum_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace, vectors.cor = vectors.cor, trace.cor = trace.cor))
    }else {
      note <- "No correction was applied to the negative eigenvalues"
      bs <- bstick.def(k3)
      bs <- c(bs, rep(0, (k - k3)))
      cum.bs <- cumsum(bs)
      res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Relative_eig", 
                         "Rel_corr_eig", "Broken_stick", "Cum_corr_eig", 
                         "Cumul_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace))
    }
  }
  class(out) <- "pcoa"
  out
}


########  function PACo ####


prepare_paco_data_test <- function (dist_A, dist_B, network) {
  
  HP <- t(as.matrix(network))
  H <- dist_A
  P <- dist_B
  # now tree_A on rows and tree_B on columns
  
  if (NROW(H) != NCOL(H)) 
    stop("dist_A should be a square matrix")
  if (NROW(P) != NCOL(P)) 
    stop("dist_B should be a square matrix")
  if (NROW(H) != NROW(HP)) {
    warning("The network matrix should have clade_A in columns. It has been translated.")
    HP <- t(HP)
  }
  if (!(NROW(H) %in% dim(HP))) 
    stop("The number of species in dist_A and network don't match")
  if (!(NROW(P) %in% dim(HP))) 
    stop("The number of species in dist_B and network don't match")
  if (!all(rownames(HP) %in% rownames(H))) 
    stop("The species names dist_A and network don't match")
  if (!all(colnames(HP) %in% rownames(P))) 
    stop("The species names dist_B and network don't match")
  H <- H[rownames(HP), rownames(HP)]
  P <- P[colnames(HP), colnames(HP)]
  HP[HP > 0] <- 1
  D <- list(H = H, P = P, HP = HP)
  class(D) <- "paco"
  return(D)
}


add_pcoord_test <- function (D) {
  HP_bin <- which(D$HP > 0, arr.ind = TRUE)
  H_PCo <- coordpcoa_test(D$H, correction = "cailliez")$vectors
  P_PCo <- coordpcoa_test(D$P, correction = "cailliez")$vectors
  D$H_PCo <- H_PCo[HP_bin[, 1,drop=F], ,drop=F]
  D$P_PCo <- P_PCo[HP_bin[, 2,drop=F], ,drop=F]
  return(D)
}

coordpcoa_test <- function (D, correction = "none", rn = NULL) {
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  bstick.def <- function(n, tot.var = 1, ...) {
    res <- rev(cumsum(tot.var/n:1)/n)
    names(res) <- paste("Stick", seq(len = n), sep = "")
    return(res)
  }
  D <- as.matrix(D)
  n <- nrow(D)
  epsilon <- sqrt(.Machine$double.eps)
  if (length(rn) != 0) {
    names <- rn
  }else {
    names <- rownames(D)
  }
  CORRECTIONS <- c("none", "lingoes", "cailliez")
  correct <- pmatch(correction, CORRECTIONS)
  if (is.na(correct)) 
    stop("Invalid correction method")
  delta1 <- centre((-0.5 * D^2), n)
  trace <- sum(diag(delta1))
  D.eig <- eigen(delta1)
  D.eig$values <- as.numeric(zapsmall(vegan::eigenvals(D.eig)))
  min.eig <- min(D.eig$values)
  zero.eig <- which(abs(D.eig$values) < epsilon)
  D.eig$values[zero.eig] <- 0
  if (min.eig > -epsilon) {
    correct <- 1
    eig <- D.eig$values
    k <- length(which(eig > epsilon))
    rel.eig <- eig[1:k]/trace
    cum.eig <- cumsum(rel.eig)
    vectors <- sweep(D.eig$vectors[, 1:k,drop=F], 2, sqrt(eig[1:k,drop=F]), FUN = "*")
    bs <- bstick.def(k)
    cum.bs <- cumsum(bs)
    res <- data.frame(eig[1:k], rel.eig, bs, cum.eig, cum.bs)
    colnames(res) <- c("Eigenvalues", "Relative_eig", "Broken_stick", 
                       "Cumul_eig", "Cumul_br_stick")
    rownames(res) <- 1:nrow(res)
    rownames(vectors) <- names
    colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                  prefix = "Axis.")
    note <- paste("There were no negative eigenvalues. No correction was applied")
    out <- (list(correction = c(correction, correct), note = note, 
                 values = res, vectors = vectors, trace = trace))
  } else {
    k <- n
    eig <- D.eig$values
    rel.eig <- eig/trace
    rel.eig.cor <- (eig - min.eig)/(trace - (n - 1) * min.eig)
    rel.eig.cor = c(rel.eig.cor[1:(zero.eig[1] - 1)], rel.eig.cor[(zero.eig[1] + 1):n], 0)
    cum.eig.cor <- cumsum(rel.eig.cor)
    k2 <- length(which(eig > epsilon))
    k3 <- length(which(rel.eig.cor > epsilon))
    vectors <- sweep(D.eig$vectors[, 1:k2,drop=F], 2, sqrt(eig[1:k2,drop=F]), 
                     FUN = "*")
    if ((correct == 2) | (correct == 3)) {
      if (correct == 2) {
        c1 <- -min.eig
        note <- paste("Lingoes correction applied to negative eigenvalues: D' = -0.5*D^2 -", 
                      c1, ", except diagonal elements")
        D <- -0.5 * (D^2 + 2 * c1)
      }else if (correct == 3) {
        delta2 <- centre((-0.5 * D), n)
        upper <- cbind(matrix(0, n, n), 2 * delta1)
        lower <- cbind(-diag(n), -4 * delta2)
        sp.matrix <- rbind(upper, lower)
        c2 <- max(Re(eigen(sp.matrix, symmetric = FALSE, 
                           only.values = TRUE)$values))
        note <- paste("Cailliez correction applied to negative eigenvalues: D' = -0.5*(D +", 
                      c2, ")^2, except diagonal elements")
        D <- -0.5 * (D + c2)^2
      }
      diag(D) <- 0
      mat.cor <- centre(D, n)
      toto.cor <- eigen(mat.cor)
      toto.cor$values <- as.numeric(zapsmall(vegan::eigenvals(toto.cor)))
      trace.cor <- sum(diag(mat.cor))
      min.eig.cor <- min(toto.cor$values)
      zero.eig.cor <- which((toto.cor$values < epsilon) & 
                              (toto.cor$values > -epsilon))
      toto.cor$values[zero.eig.cor] <- 0
      if (min.eig.cor > -epsilon) {
        eig.cor <- toto.cor$values
        rel.eig.cor <- eig.cor[1:k]/trace.cor
        cum.eig.cor <- cumsum(rel.eig.cor)
        k2 <- length(which(eig.cor > epsilon))
        vectors.cor <- sweep(toto.cor$vectors[, 1:k2,drop=F], 
                             2, sqrt(eig.cor[1:k2,drop=F]), FUN = "*")
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
        bs <- bstick.def(k2)
        bs <- c(bs, rep(0, (k - k2)))
        cum.bs <- cumsum(bs)
      }else {
        if (correct == 2) 
          cat("Problem! Negative eigenvalues are still present after Lingoes", 
              "\n")
        if (correct == 3) 
          cat("Problem! Negative eigenvalues are still present after Cailliez", 
              "\n")
        rel.eig.cor <- cum.eig.cor <- bs <- cum.bs <- rep(NA, 
                                                          n)
        vectors.cor <- matrix(NA, n, 2)
        rownames(vectors.cor) <- names
        colnames(vectors.cor) <- colnames(vectors.cor, 
                                          do.NULL = FALSE, prefix = "Axis.")
      }
      res <- data.frame(eig[1:k], eig.cor[1:k], rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Corr_eig", "Rel_corr_eig", 
                         "Broken_stick", "Cum_corr_eig", "Cum_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace, vectors.cor = vectors.cor, trace.cor = trace.cor))
    }else {
      note <- "No correction was applied to the negative eigenvalues"
      bs <- bstick.def(k3)
      bs <- c(bs, rep(0, (k - k3)))
      cum.bs <- cumsum(bs)
      res <- data.frame(eig[1:k], rel.eig, rel.eig.cor, 
                        bs, cum.eig.cor, cum.bs)
      colnames(res) <- c("Eigenvalues", "Relative_eig", 
                         "Rel_corr_eig", "Broken_stick", "Cum_corr_eig", 
                         "Cumul_br_stick")
      rownames(res) <- 1:nrow(res)
      rownames(vectors) <- names
      colnames(vectors) <- colnames(vectors, do.NULL = FALSE, 
                                    prefix = "Axis.")
      out <- (list(correction = c(correction, correct), 
                   note = note, values = res, vectors = vectors, 
                   trace = trace))
    }
  }
  class(out) <- "pcoa"
  out
}


PACo_test <- function (D, nperm = 1000, seed = NA, nullmodel = 2, symmetric = TRUE, table_biogeo = NULL) {
  
  if (!nullmodel %in% c(1,2,3)) {
    stop("Pick a null model among \"1\" or \"2\" or \"3\".")
  }
  if (!("H_PCo" %in% names(D))) {
    D <- add_pcoord_test(D)
  }
  
  suppressWarnings(
    proc <- vegan::procrustes(X = D$H_PCo, Y = D$P_PCo, symmetric = symmetric)
  )
  m2ss <- proc$ss   
  if (!is.na(seed)) set.seed(seed)
  
  
  if (nullmodel == 3) { 
    if (!all(table_biogeo$species %in% rownames(D$H_PCo))){
      stop("Provide a table_biogeo with all the species from tree_A.")
    }
    table_biogeo <- table_biogeo[match(rownames(D$H_PCo), table_biogeo$species), ]
  }
  
  ss_perm <- numeric(nperm)
  if (nullmodel == 1) {
    null_model <- vegan::nullmodel(D$HP, "r0")
    randomised_matrices <- stats::simulate(null_model, nsim = nperm)
  }
  
  for (n in 1:nperm) {
    
    if (nullmodel == 1) {
      permuted_HP <- randomised_matrices[, , n]
      permuted_HP <- permuted_HP[rownames(D$HP), colnames(D$HP)]
    }
    
    if (nullmodel == 2) { # null model 2 : random permutations of rows
      permuted_HP <- D$HP[sample(1:nrow(D$HP)), ]
      rownames(permuted_HP) <- rownames(D$HP)
    }
    
    if (nullmodel == 3) { # null model 3 : permutations constraining the biogeography
      permuted_HP <- D$HP
      for (zone in unique(table_biogeo$biogeo)) {
        w <- which(table_biogeo$biogeo==zone)
        s <- as.numeric(sample(as.character(w))) 
        permuted_HP[w,] <- as.matrix(permuted_HP[s,])
      }
      rownames(permuted_HP) <- rownames(D$HP)
    }
    
    perm_D <- list(H = D$H, P = D$P, HP = permuted_HP)
    perm_paco <- add_pcoord_test(perm_D)
    
    suppressWarnings(
      perm_proc <- vegan::procrustes(X = perm_paco$H_PCo, Y = perm_paco$P_PCo, symmetric = symmetric)
    )
    
    ss_perm[n] <- perm_proc$ss
  }
  ss_perm <- ss_perm[!is.na(ss_perm)]
  
  R2_obs <- 1 - m2ss
  R2_perm <- 1 - ss_perm
  
  pvalue <- mean(R2_perm >= R2_obs)
  
  out <- c()
  out$proc <- proc
  out$R2_obs <- R2_obs
  out$p <- pvalue
  out$R2_perm <- R2_perm
  out$ss <- m2ss
  out$nperm <- nperm
  out$nullmodel <- nullmodel
  return(out)
}


plot_H0 <- function(res, res_biogeo, method){
  
  if (method=="paco"){
    x_obs <- res$R2_obs
    x_H0 <- res$R2_perm
    x_H0_biogeo <- res_biogeo$R2_perm
    name <- "PACo R²"
    obs <- "Obs. R²"
  }
  if (method=="parafit"){
    x_obs <- res$ParaFitGlobal
    x_H0 <- res$ParaFitGlobal_perm
    x_H0_biogeo <- res_biogeo$ParaFitGlobal_perm
    name <- "ParaFit stat."
    obs <- "Obs. stat."
  }
  if (method=="mantel"){
    x_obs <- res$mantel_cor_A
    x_H0 <- res$mantel_cor_A_perm
    x_H0_biogeo <- res_biogeo$mantel_cor_A_perm
    name <- "Mantel R"
    obs <- "Obs. R"
  }
  
  x_all <- c(x_obs, x_H0, x_H0_biogeo)
  
  hist(x_H0,
       col =  adjustcolor("blue", alpha.f = 0.5),
       breaks = 20,
       freq = FALSE,
       xlim = range(x_all) + c(-0.02, 0.02) * diff(range(x_all)),
       main = "",
       xlab = name)
  
  hist(x_H0_biogeo,
       col =  adjustcolor("darkgreen", alpha.f = 0.5),
       breaks = 20,
       freq = FALSE,
       add = TRUE)
  
  abline(v = x_obs, col = adjustcolor("red", alpha.f = 0.5), lwd = 3, lty = 3)
  
  legend("topright",
         legend = c("H0", "H0 biogeo", obs),
         col = c(adjustcolor("blue", alpha.f = 0.5), adjustcolor("darkgreen", alpha.f = 0.5), adjustcolor("red", alpha.f = 0.5)),
         lwd = c(2, 2, 3),
         lty = c(1, 1, 3),
         bty = "n")
}






