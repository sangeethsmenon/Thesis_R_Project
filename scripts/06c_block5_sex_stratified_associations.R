packages <- c("tidyverse","here","splines")
need <- setdiff(packages, rownames(installed.packages()))
if (length(need)) install.packages(need, repos="https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))
suppressPackageStartupMessages({ library(jsonlite) })

if (!exists("config")) {
  cfg_rds <- here::here("results","config.rds"); cfg_json <- here::here("results","config.json")
  if (file.exists(cfg_rds)) config <- readRDS(cfg_rds) else if (file.exists(cfg_json)) config <- jsonlite::read_json(cfg_json, simplifyVector=TRUE) else stop("Config not found.")
}

b3_path <- file.path(config$paths$data_processed, "block3_scaled_pca_outliers.rds")
stopifnot(file.exists(b3_path))
b3 <- readRDS(b3_path)

dir.create(config$paths$results, recursive = TRUE, showWarnings = FALSE)

id_cols <- c("sid","sex","age_in0")
covars <- b3$covars_clean; covars$sex <- factor(covars$sex)
lev <- levels(covars$sex); stopifnot(length(lev)>=2)
idxF <- which(covars$sex == lev[1]); idxM <- which(covars$sex == lev[2])

# --- Coerce to numeric matrices (prevents 'solve' type errors) ---
to_mat <- function(df){
  M <- as.matrix(df)
  storage.mode(M) <- "double"
  M
}
metab <- to_mat(b3$metabolome_clean %>% select(-all_of(id_cols)))
prot  <- to_mat(b3$proteome_clean   %>% select(-all_of(id_cols)))

scale_cols <- function(M){
  mu <- colMeans(M)
  sdv <- apply(M,2,sd); sdv[sdv==0] <- 1
  sweep(sweep(M,2,mu,"-"),2,sdv,"/")
}

alpha <- tryCatch(config$thresholds$assoc_alpha, error=function(e) 0.05); if (is.null(alpha)) alpha <- 0.05
min_abs_r <- tryCatch(config$thresholds$assoc_min_abs_r, error=function(e) 0.10); if (is.null(min_abs_r)) min_abs_r <- 0.10

# --- Robust residualizer: QR least squares (no solve(X'X, *)) ---
make_resid_age <- function(Y, age){
  age <- as.numeric(age)
  X <- model.matrix(~ splines::ns(age, df = 3))
  # QR-based coefficients: works even if X'X is near-singular
  beta <- qr.coef(qr(X), Y)
  E <- Y - X %*% beta
  scale_cols(E)
}

M_F <- make_resid_age(metab[idxF, , drop=FALSE], covars$age_in0[idxF])
P_F <- make_resid_age(prot [idxF, , drop=FALSE], covars$age_in0[idxF])
M_M <- make_resid_age(metab[idxM, , drop=FALSE], covars$age_in0[idxM])
P_M <- make_resid_age(prot [idxM, , drop=FALSE], covars$age_in0[idxM])

cor_block <- function(P, M){
  n <- nrow(M); df <- n - 2
  Rp <- crossprod(P, M) / (n - 1); Rp[Rp>1] <- 1; Rp[Rp< -1] <- -1
  Tp <- Rp * sqrt(df / pmax(1e-12, 1 - Rp^2))
  Pp <- 2 * pt(abs(Tp), df = df, lower.tail = FALSE)
  Qp <- matrix(p.adjust(as.vector(Pp), "BH"), nrow = nrow(Pp), ncol = ncol(Pp))
  list(R=Rp, p=Pp, q=Qp, n=n, df=df)
}

resF <- cor_block(P_F, M_F)
resM <- cor_block(P_M, M_M)

prot_names  <- colnames(prot)
metab_names <- colnames(metab)
edgesF <- tibble(protein=rep(prot_names, times=length(metab_names)),
                 metabolite=rep(metab_names, each=length(prot_names)),
                 r_F=as.vector(resF$R), p_F=as.vector(resF$p), q_F=as.vector(resF$q))
edgesM <- tibble(protein=rep(prot_names, times=length(metab_names)),
                 metabolite=rep(metab_names, each=length(prot_names)),
                 r_M=as.vector(resM$R), p_M=as.vector(resM$p), q_M=as.vector(resM$q))

edges <- edgesF %>% inner_join(edgesM, by=c("protein","metabolite")) %>%
  mutate(abs_r_F=abs(r_F), abs_r_M=abs(r_M))

# Fisher z for heterogeneity (Male vs Female)
zF <- atanh(pmax(pmin(edges$r_F, 0.999999), -0.999999))
zM <- atanh(pmax(pmin(edges$r_M, 0.999999), -0.999999))
seF <- sqrt(1/(resF$n - 3))
seM <- sqrt(1/(resM$n - 3))
zdiff <- (zM - zF) / sqrt(seF^2 + seM^2)
p_diff <- 2 * pnorm(abs(zdiff), lower.tail = FALSE)
q_diff <- p.adjust(p_diff, "BH")

edges <- edges %>% mutate(p_diff = p_diff, q_diff = q_diff, r_diff = r_M - r_F)

edges_F_sig    <- edges %>% filter(q_F <= alpha, abs_r_F >= min_abs_r)
edges_M_sig    <- edges %>% filter(q_M <= alpha, abs_r_M >= min_abs_r)
edges_diff_sig <- edges %>% filter(q_diff <= alpha, abs(r_diff) >= 0.10)

summ <- tibble(
  view = "sex_stratified_age_adjusted",
  n_F = resF$n, n_M = resM$n,
  n_tests = nrow(edges),
  alpha = alpha, min_abs_r = min_abs_r,
  F_sig = nrow(edges_F_sig),
  M_sig = nrow(edges_M_sig),
  diff_sig = nrow(edges_diff_sig)
)

readr::write_csv(summ, file.path(config$paths$results, "block5c_sex_stratified_summary.csv"))
readr::write_csv(edges_F_sig %>% arrange(q_F, desc(abs_r_F)),
                 file.path(config$paths$results, "block5c_edges_F_q_le_alpha_absr.csv"))
readr::write_csv(edges_M_sig %>% arrange(q_M, desc(abs_r_M)),
                 file.path(config$paths$results, "block5c_edges_M_q_le_alpha_absr.csv"))
readr::write_csv(edges %>% arrange(q_diff) %>% slice_head(n=20000),
                 file.path(config$paths$results, "block5c_edges_top20k_by_pdiff.csv"))

saveRDS(list(edges=edges, summary=summ,
             n_F=resF$n, n_M=resM$n),
        file = file.path(config$paths$data_processed, "block5c_sex_stratified_associations.rds"))

cat("Block 5C (sex-stratified, age-adjusted) complete.\n")
