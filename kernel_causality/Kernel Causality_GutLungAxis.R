rm(list = ls())

if(!"generalCorr"%in%installed.packages()){
  install.packages("generalCorr", repos="http://cran.case.edu/")
}
suppressPackageStartupMessages({
  library(generalCorr)
  library(dplyr); library(readr); library(stringr); library(scales); library(tidyr)
})

library(readxl)

if(!"mutoss"%in%installed.packages()){
  install.packages("mutoss", repos="https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(mutoss))

#####################################################################################

## 2) 파일 지정 및 전처리

setwd('C:/Users/tr432/Downloads/Metabolomics_GutLungAxis_re-pipeline/merged_correlation')

df <- read_excel('merged_corr_input.xlsx', sheet='merged_kernel')#,
#                 show_col_types = FALSE)


## (선택1) 특정 변수만 사용
## 1) kernel에 넣을 변수 이름
vars_use <- c(
  'g__Bacteroides_H_857956', 
  'g__Bifidobacterium_388775', 
  'g__Escherichia', 
  'g__Lactobacillus', 
  'g__Phocaeicola_A', 
  'g__Prevotella', 
  'g__Proteus', 
  'Lactobacillus_jensenii', 
  'Bacillus_subtilis', 
  'Staphylococcus_aureus', 
  'Phocaeicola_vulgatus', 
  'Bacteroides_thetaiotaomicron', 
  'Enterobacter_cancerogenus', 
  'Enterobacter_cloacae_complex_sp.', 
  'Enterobacter_soli', 
  'Enterobacter_mori', 
  'Yersinia_pestis', 
  'Morganella_morganii', 
  'BM_whole', 
  'BM_wbc', 
  'Spleen', 
  'b.w', 
  'Spleen/b.w', 
  'Fecal_205.09707-5.741', 
  'Fecal_192.06554-7.209', 
  'Fecal_190.0499-5.9', 
  'Fecal_176.07063-8.668', 
  'Fecal_161.10737-6.299', 
  'Fecal_224.0553-2.666', 
  'Fecal_221.09201-2.066', 
  'Fecal_192.0655-13.16', 
  'Fecal_177.10222-6.754', 
  'Fecal_206.04465-6.374', 
  'Fecal_219.11272-7.086', 
  'Fecal_118.06507-10.451', 
  'Fecal_160.07568-5.773', 
  'Fecal_61.02846-1.384', 
  'Fecal_75.04402-1.443', 
  'Fecal_150.05836-1.457', 
  'Fecal_206.08106-8.623', 
  'Fecal_190.08619-9.491', 
  'Fecal_104.10688-1.024', 
  'Fecal_496.34005-14.377', 
  'Fecal_522.35578-14.691', 
  'Fecal_146.11762-1.098', 
  'Fecal_508.37662-16.157', 
  'Fecal_518.3247-13.08', 
  'Fecal_466.29327-14.068', 
  'Fecal_780.5541-17.92', 
  'Serum_188.07062-5.731', 
  'Serum_176.07059-8.684', 
  'Serum_162.09137-7.84', 
  'Serum_215.01613-1.47', 
  'Serum_61.02844-1.423', 
  'Serum_150.05836-1.455', 
  'Serum_190.08624-9.5', 
  'Serum_496.33995-14.389', 
  'Serum_522.35575-14.699', 
  'Serum_104.10687-1.017', 
  'Serum_526.29348-13.683', 
  'Serum_542.32483-13.118', 
  'Serum_518.32462-13.233'
  )

## (선택2) 모든 변수 사용
vars_use <- df %>%
  dplyr::pull(Features) %>%
  as.character() %>%
  unique()

## 2) Features에서 필요한 것만 행으로 선택
dat0 <- df %>%
  dplyr::filter(Features %in% vars_use)

## 3) Features를 행이름으로, 나머지는 샘플 값만 남김
mat <- dat0 %>%
  dplyr::select(-Features) %>%
  as.data.frame()
rownames(mat) <- dat0$Features   # 행이름 = 변수명

## 4) 전치: 행 = 샘플, 열 = 변수
dat <- as.data.frame(t(mat))
colnames(dat) <- rownames(mat)   # 열 이름 = 변수명들(g_Phocaeciola_A, ...)

## 5) 숫자형 변환만 수행
dat_num <- dat %>%
  dplyr::mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

## pairwise complete 방식으로 Kernel causality 수행
options(np.messages = FALSE)

vars <- colnames(dat_num)

MIN_N <- 2   # 변수쌍별 최소 샘플 수. 샘플 넉넉하면 6 추천

res_list <- list()
k <- 1

for (i in 1:(length(vars) - 1)) {
  for (j in (i + 1):length(vars)) {
    
    xname <- vars[i]
    yname <- vars[j]
    
    ## 해당 변수쌍에서만 NA 제거
    tmp <- dat_num[, c(xname, yname), drop = FALSE] %>%
      tidyr::drop_na()
    
    ## 둘 다 값 있는 샘플 수가 너무 적으면 제외
    if (nrow(tmp) < MIN_N) next
    
    ## 값이 전부 동일한 변수는 correlation/kernel 계산 불가
    if (sd(tmp[[1]], na.rm = TRUE) == 0 || sd(tmp[[2]], na.rm = TRUE) == 0) next
    
    ## Kernel Causality (generalCorr) -> allPairs 실행
    out <- tryCatch({
      suppressWarnings(
        as.data.frame(allPairs(tmp, typ = 3, dig = 4))
      )
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(out) && nrow(out) > 0) {
      out$n_used <- nrow(tmp)
      out$X_input <- xname
      out$Y_input <- yname
      
      res_list[[k]] <- out
      k <- k + 1
    }
  }
}

m1_df <- dplyr::bind_rows(res_list)

#####################################################################################

## p-value 컬럼 찾기
pcol <- intersect(c("p-val","p_val","p.value","pvalue","P"), names(m1_df))
if(length(pcol) == 0){
  stop("allPairs 결과에서 p-value 컬럼을 못 찾았습니다. 현재 컬럼: ",
       paste(names(m1_df), collapse=", "))
}
pcol <- pcol[1]


## BKY(q) 계산 -------------------------------------------------------------

## 빠른 BH adjusted p-value 함수
bh_adjust_fast <- function(p) {
  p <- as.numeric(p)
  q <- rep(NA_real_, length(p))
  
  keep <- !is.na(p) & is.finite(p)
  p0 <- p[keep]
  
  m <- length(p0)
  if (m == 0) return(q)
  
  o <- order(p0)
  ro <- order(o)
  
  p_sorted <- p0[o]
  q_sorted <- p_sorted * m / seq_len(m)
  q_sorted <- rev(cummin(rev(q_sorted)))
  q_sorted <- pmin(q_sorted, 1)
  
  q[keep] <- q_sorted[ro]
  q
}

## 빠른 BKY two-stage FDR 함수
bky_adjust_fast <- function(p, alpha = 0.05) {
  p <- as.numeric(p)
  q <- rep(NA_real_, length(p))
  
  keep <- !is.na(p) & is.finite(p)
  p0 <- p[keep]
  
  m <- length(p0)
  if (m == 0) return(q)
  
  ## 1단계: alpha' = alpha / (1 + alpha)
  alpha_prime <- alpha / (1 + alpha)
  
  q_bh <- bh_adjust_fast(p0)
  r1 <- sum(q_bh <= alpha_prime, na.rm = TRUE)
  
  if (r1 == 0 || r1 == m) {
    q0 <- q_bh * (1 + alpha)
  } else {
    m0 <- m - r1
    q0 <- q_bh * (m0 / m) * (1 + alpha)
  }
  
  q0 <- pmin(q0, 1)
  q[keep] <- q0
  
  q
}


## 적용
p_vec <- as.double(m1_df[[pcol]])
m1_df$q_bky <- bky_adjust_fast(p_vec, alpha = 0.05)


## 저장
readr::write_csv(m1_df,  "kernel_allPairs_BKY.csv")

m1_sig <- dplyr::filter(m1_df, !is.na(`p-val`) & `p-val` < 0.05)  #Edit: `p-val` or q_bky
#readr::write_csv(m1_sig, "kernel_allPairs_BKY_qlt0.05.csv")

#####################################################################################

## 4) Cytoscape용 엣지/노드 테이블 생성

edges_src <- m1_sig  # 또는 read_tsv(m1_path, show_col_types = FALSE)

abs_r_cut <- 0.0          # 약한 엣지 제거 임계값
width_range   <- c(1, 5)       # Cytoscape Web에서 보이는 엣지 두께(px)
top_n_nodes   <- 999            # 상위 노드 개수

edges_full <- edges_src %>%
  transmute(
    X, Y, Cause,
    r    = as.double(r),
    r_xy = as.double(`r*x|y`),
    r_yx = as.double(`r*y|x`)
  ) %>%
  mutate(
    abs_r = abs(r),
    Source    = if_else(Cause == X, X, Y),
    Target    = if_else(Cause == X, Y, X),
    Weight    = abs_r
  ) %>%
  filter(Source != Target, abs_r >= abs_r_cut)

node_strength <- bind_rows(
  transmute(edges_full, Node = Source, w = Weight),
  transmute(edges_full, Node = Target, w = Weight)
) %>%
  group_by(Node) %>%
  summarise(strength = sum(w, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(strength))

top_nodes <- head(node_strength$Node, top_n_nodes)

edges_top <- edges_full %>%
  filter(Source %in% top_nodes, Target %in% top_nodes) %>%
  mutate(
    Interaction = "causes",
    Polarity    = if_else(r >= 0, "pos", "neg"),
    EdgeWidth   = round(rescale(Weight, to = width_range), 2)
  ) %>%
  dplyr::select(Source, Target, Interaction, EdgeWidth, Polarity)

nodes_top <- tibble::tibble(Node = sort(unique(c(edges_top$Source, edges_top$Target))))

# 저장
readr::write_tsv(edges_top, "edges_kernel.tsv")
readr::write_tsv(nodes_top, "nodes_kernel.tsv")

#####################################################################################

### Cytoscape 스타일 메모

# Edges
#  Label Font Size : 12
#  Line Type : Column-Polarity / Mapping Type-discrete
#              pos-solid / neg-dashed
#  Stroke Color : #C9C9C9 (Edge color to arrows)
#  Target Arrow Shape : Column-Interaction / Mapping Type-discrete
#                       causes-triangle
#  Width : Column-EdgeWidth / Mapping Type-passthrough


# Nodes
#  Border Color : #C9C9C9
#  Border Width : 3
#  Fill Color : #FFB000
#  Height : 13 (Lock node width and height)
#  Label : Column-name / Mapping Type-passthrough
#  Label Font Size : 15
#  Label Position : bottom-right
#  Opacity : 80%

#####################################################################################

