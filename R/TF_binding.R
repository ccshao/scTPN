#' Identify mutual bindings
#'
#' Examine the mutual bindings for each pair of TF in each other's promoter sequences.
#'
#' @param mtx a symmetric TFs connectivity matrix.
#' @param species in which species looking for the pfm of TFs, e.g. "Hsapiens" (default), "Mmusculus".
#' @param gene_db org.Xx.eg.db annotation converting gene symbol to entrez id, e.g., \code{org.Mm.eg.db::org.Mm.egSYMBOL2EG}.
#' @param bsgenome bsgenome of the studied species, e.g., "BSgenome.Hsapiens.UCSC.hg19".
#' @param txdb txdb database, e.g., "TxDb.Hsapiens.UCSC.hg19.knownGene" for human.
#' @param upstream upstream range of the TSS.
#' @param downstream downstream range of the TSS.
#' @param min_score min_score for the PWM searching used in \code{matchPWM} (default "85\%").
#' @param verbose whether show more info.
#'
#' @details
#' The input is the TFs netowrk , predicted from \code{skel_net} or \code{pcn}.
#' The binding sites of a pair of TFs were predicted on both forward
#' and reverse strand for available PWMs sepecfied by parameters in a mutual way (i.e., pwms of A on the promter seqs of B, and the other way around).
#' The numbers are TF PWMs bound on promoter regions (one TF might have multiple PWMs), but not the number of actual binding sites.
#'
#' The TF position weight matrices were obainted from MotifDb, and binding sites are identified by \code{matchPWM}.
#' @return A matrix indicates whether there are TFs binding sites between genes (0 or number of bounded PWMs).
#' @export
pwm_pred <- function(mtx,
    species = "Hsapiens",
    gene_db,
    bsgenome,
    txdb,
    upstream = 1000,
    downstream = 200,
    min_score = "85%",
    verbose = FALSE) {
  rn  <- i <- j <- NULL

  if (!is.matrix(mtx)) stop("(EE) Input should be a matrix.")

  #- Only take the upper tri as mtx is symm
  mtx[lower.tri(mtx)] <- 0
  inds   <- which(mtx %!=% 0, arr.ind = TRUE)
  idx_dt <- data.table(rn = rownames(mtx)[inds[, 1]], yn = colnames(mtx)[inds[, 2]])[order(rn)]

  #- An empty matrix.
  tf_mtx   <- mtx
  tf_mtx[] <- 0

  #- Get seqs.
  all_tf    <- unique(c(idx_dt$rn, idx_dt$yn))
  all_trans <- GenomicFeatures::transcriptsBy(txdb, by = "gene")
  tf_seq    <- sapply(all_tf, fn_seq, gene_db, bsgenome, all_trans, upstream, downstream, simplify = FALSE)
  tf_seq[lengths(tf_seq) == 0] <- NULL

  message(sprintf("(II) Promoter seqs found for %i/%i TFs.", length(tf_seq), nrow(mtx)))

  if (length(tf_seq) == 0) stop("(EE) NO promter seq found for TFs! Check whether ids are correct.")

  #- Get pwms. Note one TFs might have multiple pwms.
  tf_pfm <- sapply(all_tf, function(x, species) suppressMessages(as.list(MotifDb::query(MotifDb::MotifDb, c(x, species)))), species, simplify = FALSE)
  tf_pfm[lengths(tf_pfm) == 0] <- NULL

  message(sprintf("(II) PWMs found for %i/%i TFs.", length(tf_pfm), nrow(mtx)))

  if (length(tf_pfm) == 0) stop("(EE) NO pwm found for TFs! Check whether MotifDb has those TFs.")

  res <- foreach(i = idx_dt$rn, j = idx_dt$yn) %do% {
    fn_tf_pair_search(i, j, tf_pfm, tf_seq, min_score, verbose = verbose)
  } %>%
  as.data.table %>%
  t

  tf_mtx[res[, 1:2]] <- as.numeric(res[, 3])

  #- Make symmetric matrix
  tf_mtx[lower.tri(tf_mtx)] <- t(tf_mtx)[lower.tri(tf_mtx)]
  return(tf_mtx)
}

fn_seq <- function(x, gene_db, bsgenome, all_trans, upstream, downstream) {
  symbol <- gene_id <- NULL

  target_id <- AnnotationDbi::mget(x, gene_db, ifnotfound = NA)[[1]]

  #- Take the promoter of first transcript.
  if (!is.na(target_id) && (target_id %in% names(all_trans))) {
    promseq <- GenomicFeatures::getPromoterSeq(all_trans[[target_id]][1, ], bsgenome, upstream = upstream, downstream = downstream) %>% unlist
  }
}

fn_tf_pair_search <- function(tf_1, tf_2, tf_pfm, tf_seq, min_score, verbose = FALSE) {
  #- Total binding pfw found on the promoter regions with swapping targets.
  res_1 <- fn_search(tf_1, tf_2, tf_pfm, tf_seq, min_score, verbose = verbose)
  res_2 <- fn_search(tf_2, tf_1, tf_pfm, tf_seq, min_score, verbose = verbose)
  #- Sum of TF pfm having hits in either of directions.
  return(c(tf_1, tf_2, res_1 + res_2))
}

#- Here
fn_search <- function(in_tf, target_tf, tf_pfm, tf_seq, min_score, verbose = FALSE) {
  #- Return 0 if no TF pfm or no hits, or return how many pfm could be identified summaize across TFs,
  #- but not the actuall binding sites.
  n_pwm <- 0

  if (in_tf %in% names(tf_pfm) && (target_tf %in% names(tf_seq))) {
    #- Identify hits for each pfm.
    if (verbose) message("(II) searching pwms of ", in_tf, " in the promoter of ", target_tf, ".")
    res <- vapply(tf_pfm[[in_tf]], pwm_match, numeric(1), tf_seq[[target_tf]], min_score)
    #- Record how many pwms found hits in promoter seq.
    n_pwm <- sum(res > 0)
  } else {
    message("(II) **NO** search. ", in_tf,  " pwm info: ", in_tf %in% names(tf_pfm), "; ", target_tf, " seq info: ", target_tf %in% names(tf_seq), ".")
  }
  return(n_pwm)
}

#' Get TF binding sites
#'
#' Predict the binding sites of a TF on a target sequence by searching pwm.
#'
#' @param tf_pfm position frequency matrix, obtained from MotifDb.
#' @param target_seq target gene promoter seq, a DNSstringset object.
#' @param min_score min_score for the PWM searching used in \code{matchPWM} (default "85\%").
#' @return The number of pfm binding sites on both postive and negative strands.
#' @export
pwm_match <- function(tf_pfm, target_seq, min_score = "85%") {
  tf_pwm  <- round(100 * tf_pfm)
  #- Forward and reverse strand bindings.
  n_b_1  <- Biostrings::matchPWM(tf_pwm, target_seq, min.score = min_score)
  n_b_2  <- Biostrings::matchPWM(Biostrings::reverseComplement(tf_pwm), target_seq, min.score = min_score)
  n_bdi <- length(n_b_1) + length(n_b_2)
  return(n_bdi)
}
