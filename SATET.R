#Source Saddlepoint Approximation functions
source("~/SPA_functions.R")

################ MAIN FUNCTION ###################

SATET <- function(Gmat_case,Gmat_ctrl,
                  domain_snp_dt,
                  agg_snps=FALSE,
                  thresh_quant="ID",
                  thresh_val,
                  glm_input=NULL,
                  teststat="FET",
                  midp=TRUE,
                  score.cutoff=2,
                  use.marg.p=TRUE,
                  use.SPA.score=TRUE,
                  L=3,
                  alpha=0.05){
  
  require(purrr)
  require(data.table)
  require(dplyr)
  
  #Input Parameters
  in.params <- list(agg_snps=FALSE,
                    thresh_quant="ID",
                    thresh_val=thresh_val,
                    glm_input=NULL,
                    teststat="FET",
                    midp=TRUE,
                    score.cutoff=2,
                    use.marg.p=TRUE,
                    use.SPA.score=TRUE,
                    L=L,
                    alpha=alpha)
  
  #Global objects
  p.hats <- rep(NA,L)
  S.l.list = vector("list",L)
  S.list=vector("list",L)
  S.list.dc=vector("list",L)
  LayerSumm = vector("list",L)
  
  FD_approx_prev = 0
  D_approx_prev = 0
  
  #Some processing
  if(class(Gmat_case)%in%c("CsparseMatrix","ngCMatrix")){
    mat1 <- unname(as(Gmat_case,"dgCMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgCMatrix"))
  } else{
    mat1 <- unname(as(Gmat_case,"dgTMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgTMatrix"))
  }
  
  mat_all <- rbind(mat1,mat0)
  
  N1 <- nrow(mat1)
  N0 <- nrow(mat0)
  N <- N1 + N0
  
  struct_map <- domain_snp_dt
  
  if(agg_snps){### With Aggregation
    
    if(thresh_quant=="ID"){
      struct_map <- struct_map %>%
        mutate_at("ID",as.character) %>% #Coerce ID from factor to character
        group_by(domainID) %>%
        mutate(lastObsFlagDomain = as.integer(row_number() == n())) %>%
        mutate(num_unique = compute_counts(thresh_val, ID, lastObsFlagDomain)) %>%
        ungroup() %>%
        mutate(group = cumsum(c(-1L, diff(num_unique)) <= 0L)) %>%
        group_by(domainID) %>%
        mutate(group2 = ifelse(group==max(group) & last(num_unique) < thresh_val,
                               max(group)-1L,group)) %>%
        ungroup() %>%
        mutate(L1 = rleid(group2)) %>%
        dplyr::select(-c(contains("group"),lastObsFlagDomain,num_unique)) %>%
        mutate_at("ID",as.factor) %>% #coerce ID back to factor
        data.table()
      
    } else if(thresh_quant=="MAC"){
      struct_map <- struct_map %>%
        group_by(domainID) %>%
        mutate(cumsum_val = sum_reset_at(mac, thresh_val)) %>%
        mutate(next_group1 = ifelse(lag(cumsum_val) >= thresh_val | row_number() == 1, 1, 0)) %>% ## binary interpretation of whether there should be a new group
        ungroup %>%
        mutate(group1 = cumsum(next_group1)) %>%
        ungroup %>%
        group_by(domainID) %>%
        mutate(group2 = ifelse(group1==max(group1) & last(mac) < thresh_val,
                               max(group1)-1,group1)) %>%
        ungroup() %>%
        mutate(L1 = rleid(group2)) %>%
        select(-contains("group")) %>%
        data.table()
    } else{ #thresh_quant=="numvar"
      
      #struct_map$col_maf <- Matrix::colMeans(mat_all)
      struct_map$col_nvar <- Matrix::colSums(mat_all)
      struct_map$col_case_nvar <- Matrix::colSums(mat1)
      
      struct_map <- struct_map %>%
        group_by(domainID) %>%
        mutate(cumsum_val = sum_reset_at(col_nvar, thresh_val)) %>%
        mutate(next_group1 = ifelse(lag(cumsum_val) >= thresh_val | row_number() == 1, 1, 0)) %>% ## binary interpretation of whether there should be a new group
        ungroup %>%
        mutate(group1 = cumsum(next_group1)) %>%
        ungroup %>%
        group_by(domainID) %>%
        mutate(group2 = ifelse(group1==max(group1) & last(col_nvar) < thresh_val,
                               max(group1)-1,group1)) %>%
        ungroup() %>%
        mutate(L1 = rleid(group2)) %>%
        select(-contains("group")) %>%
        data.table()
    }
    
    
    #Leaf attribute matrix
    leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
    
    #Calculate test statistics
    if(teststat=="FET"){
      #For FET statistic
      
      leaf_mat_all@x <- ifelse(leaf_mat_all@x>0,1,0) #binarize
      
      case_colSums = Matrix::colSums(leaf_mat_all[1:N1,])
      all_colSums = Matrix::colSums(leaf_mat_all)
      
      pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                      case_colSums = case_colSums,
                                      all_colSums = all_colSums,
                                      midp=midp)
      
    } else{
      #For score statistic
      
      if(!is.null(glm_input) && ncol(glm_input)>1){
        #With covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=glm_input[,1],cov=glm_input[,-1],
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=glm_input[,1],cov=glm_input[,-1],minmac=1)
        }
        
      } else{
        #Without covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,minmac=1)
        }
        
      }
      
      score.1 <- score.test$Tstat.sign
      pvals.1 <- score.test$p.value
    }
    
  } else{#### No Aggregation
    #Treat each snp as a bottom-layer leaf
    struct_map <- struct_map[, L1 := snpID]
    
    if(teststat=="FET"){
      #Get marginals - qualifying variants for cases and all
      case_colSums = Matrix::colSums(mat1)
      all_colSums = Matrix::colSums(mat_all)
      
      pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                      case_colSums = case_colSums,
                                      all_colSums = all_colSums,
                                      midp=midp)
      
    } else{
      #For score statistic
      leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
      
      if(!is.null(glm_input) && ncol(glm_input)>1){
        #With covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=glm_input[,1],cov=glm_input[,-1],
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=glm_input[,1],cov=NULL,minmac=1)
        }
        
      } else{
        #Without covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,minmac=1)
        }
        
      }
      
      score.1 <- score.test$Tstat.sign
      pvals.1 <- score.test$p.value
    }
    
  }
  
  #Copy object - this is not changed
  struct_map_ <- struct_map
  
  #If NaN, set to zero
  #(occurs when genotype vector is zero, so denominator of score stat is zero)
  pvals.1[is.na(pvals.1)] <- 1
  
  rm(mat1,mat0,mat_all)
  
  for(l in seq(L)){
    
    if(l==1){
      
      m.l <- length(pvals.1)
      
      pvals.l = pvals.1
      
      p.hat.l = est.p.hat(l=l,
                          D_prev = D_approx_prev,
                          FD_prev = FD_approx_prev,
                          pvals_l=pvals.l,
                          alpha=alpha)
      
      S.l <- which(pvals.l <= p.hat.l) 
      
      if(teststat=="score"){
        LayerSumm[[l]]=list("score.l"=score.1,
                            "pvals.l"=pvals.l,
                            "S.l"=S.l,
                            "struct_map"=struct_map)
      } else{
        LayerSumm[[l]]=list("pvals.l"=pvals.l,
                            "S.l"=S.l,
                            "struct_map"=struct_map)
      }
      
    }else { #l>1
      
      #Update struct_map
      #First remove domains with fewer than 2^(l-1) L1 leaves
      struct_map <- struct_map[ , if (.N >= 2^(l-1) && uniqueN(L1) >= 2^(l-1)) .SD, by = .(domainID)]
      
      #Create layer l indices based on grouping remaining L1 indices by domain
      Ll <- paste0("L",l)
      struct_map <-  struct_map[, (Ll) := {
        ## modify rle values
        x <- ceiling(rleid(L1) / 2^(l-1))
        n <- uniqueN(L1)
        if(n > 1 && between(n %% 2^(l-1),1,2^(l-1)-1)) {
          x[x == x[.N]] <- x[.N] - 1
        }
        x
      }, by = .(domainID)][, (Ll) := rleid(domainID, get(Ll))] #reassign groups by domain and Ll
      
      
      #Layer l-1 information
      Llm1 <- paste0("L",l-1)
      struct_map_tmp <- struct_map[, .(numSNP_lev = uniqueN(snpID),
                                       numL1_lev = uniqueN(L1),
                                       numLlm1_lev = uniqueN(get(Llm1))),
                                   by = .(domainID,get(Ll))]
      setnames(struct_map_tmp,old="get",new=Ll)
      
      #Number of layer 1 leaves used to create each layer l leaf
      numL1_leaves <- struct_map_tmp$numL1_lev
      
      numLlm1_leaves <- struct_map_tmp$numLlm1_lev
      
      #Subset layer 1 pvals to reflect current leaf indices
      pvals.sub = pvals.lm1[unique(struct_map[,get(Llm1)])]
      
      #Transform layer l-1 p-values before aggregating
      Q.sub = -2*log(pvals.sub)
      
      #Aggregate L(l-1) leaves according to partition defined on the lth layer
      Q.l <- as.vector(tapply(Q.sub,rep(seq_along(numLlm1_leaves),
                                        times=numLlm1_leaves),sum))
      
      m.l = length(Q.l)
      
      #Convert Q.l's to marginal p-values
      #Might be sub-uniform
      if(use.marg.p){
        pvals.l <- vapply(seq(m.l), function(x) pchisq(Q.l[x],df=2*numLlm1_leaves[x],
                                                       lower.tail = FALSE), numeric(1))
      } else{
        #Calculate conditional p-values depending on number of leaves merged from layer l-1
        pvals.l <- vapply(seq(m.l), function(x) ifelse(numLlm1_leaves[x]==2,
                                                       condChiSqProb2(Q.l[x],t_prev=-2*log(p.hats[l-1]),
                                                                      nu_prev_vec=c(2,2)),
                                                       condChiSqProb3(Q.l[x],t_prev=-2*log(p.hats[l-1]),
                                                                      nu_prev_vec=c(2,2,2))),
                          numeric(1))
      }
      
      
      #Obtain layer specific-threshold
      p.hat.l = est.p.hat(l=l,
                          D_prev = D_approx_prev,
                          FD_prev = FD_approx_prev,
                          pvals_l=pvals.l,
                          alpha=alpha)
      
      
      S.l <- which(pvals.l <= p.hat.l) #w.r.t. layer l indices
      
      
      LayerSumm[[l]]=list("Q.l"=Q.sub,
                          "pvals.l"=pvals.l,
                          "S.l"=S.l,
                          "struct_map"=struct_map)
      
    }
    
    
    #### Post-layer processing
    
    if(l==1){
      S = S.l
    } else{
      R.l <- unique(struct_map$L1[struct_map[,get(Ll)]%in%S.l])
      S = c(S,R.l)
    }
    
    #Save pvalues from previous layer
    pvals.lm1 = pvals.l
    
    p.hats[l] = p.hat.l
    S.list[[l]] = S #discovered layer 1 leaves
    S.list.dc[[l]] = unique(struct_map_$domain[struct_map_$L1%in%S]) #vector of domain indices
    
    FD_approx_prev = FD_approx_prev + m.l*p.hats[l]
    D_approx_prev = length(S.l)
    
    #Update struct_map
    struct_map <- struct_map[!(struct_map$L1 %in% S), ]
    
  }
  
  
  return(list("params"=in.params,
              "p.hats"=p.hats,
              "S.list"=S.list,
              "S.list.dc"= S.list.dc,
              "LayerSumm"=LayerSumm))
}

#With rank aggregation
SATET_rank <- function(Gmat_case,Gmat_ctrl,
                       domain_snp_dt,
                       agg_snps=FALSE,
                       thresh_quant="ID",
                       thresh_val,
                       glm_input=NULL,
                       teststat="FET",
                       midp=TRUE,
                       score.cutoff=2,
                       use.marg.p=TRUE,
                       use.SPA.score=TRUE,
                       L=3,
                       alpha=0.05){
  
  require(purrr)
  require(data.table)
  require(dplyr)
  
  #Input Parameters
  in.params <- list(agg_snps=agg_snps,
                    thresh_quant=thresh_quant,
                    thresh_val=thresh_val,
                    glm_input=glm_input,
                    teststat=teststat,
                    midp=midp,
                    score.cutoff=score.cutoff,
                    use.marg.p=use.marg.p,
                    use.SPA.score=use.SPA.score,
                    L=L,
                    alpha=alpha)
  
  #Global objects
  p.hats <- rep(NA,L)
  S.l.list = vector("list",L)
  S.list=vector("list",L)
  S.list.dc=vector("list",L)
  LayerSumm = vector("list",L)
  
  FD_approx_prev = 0
  D_approx_prev = 0
  
  #Some processing
  if(class(Gmat_case)%in%c("CsparseMatrix","ngCMatrix")){
    mat1 <- unname(as(Gmat_case,"dgCMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgCMatrix"))
  } else{
    mat1 <- unname(as(Gmat_case,"dgTMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgTMatrix"))
  }
  
  mat_all <- rbind(mat1,mat0)
  
  N1 <- nrow(mat1)
  N0 <- nrow(mat0)
  N <- N1 + N0
  
  struct_map <- domain_snp_dt
  
  if(agg_snps){### With Aggregation
    
    if(thresh_quant=="ID"){
      struct_map <- struct_map %>%
        mutate_at("ID",as.character) %>% #Coerce ID from factor to character
        group_by(domainID) %>%
        mutate(lastObsFlagDomain = as.integer(row_number() == n())) %>%
        mutate(num_unique = compute_counts(thresh_val, ID, lastObsFlagDomain)) %>%
        ungroup() %>%
        mutate(group = cumsum(c(-1L, diff(num_unique)) <= 0L)) %>% #need <= to handle case of exactly 1 leaf in a domain
        group_by(domainID) %>%
        mutate(group2 = ifelse(group==max(group) & last(num_unique) < thresh_val,
                               max(group)-1L,group)) %>%
        ungroup() %>%
        mutate(L1 = rleid(group2)) %>%
        dplyr::select(-c(contains("group"),lastObsFlagDomain,num_unique)) %>%
        mutate_at("ID",as.factor) %>% #coerce ID back to factor
        data.table()
      
    } else if(thresh_quant=="MAC"){
      struct_map <- struct_map %>%
        group_by(domainID) %>%
        mutate(cumsum_val = sum_reset_at(mac, thresh_val)) %>%
        mutate(next_group1 = ifelse(lag(cumsum_val) >= thresh_val | row_number() == 1, 1, 0)) %>% ## binary interpretation of whether there should be a new group
        ungroup %>%
        mutate(group1 = cumsum(next_group1)) %>%
        ungroup %>%
        group_by(domainID) %>%
        mutate(group2 = ifelse(group1==max(group1) & last(mac) < thresh_val,
                               max(group1)-1,group1)) %>%
        ungroup() %>%
        mutate(L1 = rleid(group2)) %>%
        select(-contains("group")) %>%
        data.table()
    } else{ #thresh_quant=="numvar"
      
      struct_map$col_nvar <- Matrix::colSums(mat_all)
      struct_map$col_case_nvar <- Matrix::colSums(mat1)
      
      struct_map <- struct_map %>%
        group_by(domainID) %>%
        mutate(cumsum_val = sum_reset_at(col_nvar, thresh_val)) %>%
        mutate(next_group1 = ifelse(lag(cumsum_val) >= thresh_val | row_number() == 1, 1, 0)) %>% ## binary interpretation of whether there should be a new group
        ungroup %>%
        mutate(group1 = cumsum(next_group1)) %>%
        ungroup %>%
        group_by(domainID) %>%
        mutate(group2 = ifelse(group1==max(group1) & last(col_nvar) < thresh_val,
                               max(group1)-1,group1)) %>%
        ungroup() %>%
        mutate(L1 = rleid(group2)) %>%
        select(-contains("group")) %>%
        data.table()
    }
    
    #Leaf attribute matrix
    leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
    
    #Calculate test statistics
    if(teststat=="FET"){
      #For FET statistic
      
      leaf_mat_all@x <- ifelse(leaf_mat_all@x>0,1,0) #binarize
      
      case_colSums = Matrix::colSums(leaf_mat_all[1:N1,])
      all_colSums = Matrix::colSums(leaf_mat_all)
      
      pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                      case_colSums = case_colSums,
                                      all_colSums = all_colSums,
                                      midp=midp)
      
      signs.FET.1 <- ifelse(case_colSums > all_colSums-case_colSums,1,-1)
      
      signed.pvals.1 <- signs.FET.1*pvals.1
      
    } else{
      #For score statistic
      
      if(!is.null(glm_input) && ncol(glm_input)>1){
        #With covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=glm_input[,1],cov=glm_input[,-1],
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=glm_input[,1],cov=glm_input[,-1],minmac=1)
        }
        
      } else{
        #Without covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,minmac=1)
        }
        
      }
      
      score.1 <- score.test$Tstat.sign
      pvals.1 <- score.test$p.value
      
      signed.pvals.1 <- sign(score.1)*pvals.1
    }
    
  } else{#### No Aggregation
    #Treat each snp as a bottom-layer leaf
    struct_map <- struct_map[, L1 := snpID]
    
    if(teststat=="FET"){
      #Get marginals - qualifying variants for cases and all
      case_colSums = Matrix::colSums(mat1)
      all_colSums = Matrix::colSums(mat_all)
      
      pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                      case_colSums = case_colSums,
                                      all_colSums = all_colSums,
                                      midp=midp)
      
    } else{
      #For score statistic
      leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
      
      if(!is.null(glm_input) && ncol(glm_input)>1){
        #With covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=glm_input[,1],cov=glm_input[,-1],
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=glm_input[,1],cov=NULL,minmac=1)
        }
        
      } else{
        #Without covariates
        if(use.SPA.score){
          score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                                 pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,
                                                 minmac=1,Cutoff=score.cutoff)
        } else{
          score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                         pheno=rep(c(1,0),times=c(N1,N0)),cov=NULL,minmac=1)
        }
        
      }
      
      score.1 <- score.test$Tstat.sign
      pvals.1 <- score.test$p.value
      
      signed.pvals.1 <- sign(score.1)*pvals.1
    }
    
  }
  
  #Copy object - this is not changed
  struct_map_ <- struct_map
  
  #If NaN, set to zero
  #(occurs when genotype vector is zero, so denominator of score stat is zero)
  pvals.1[is.na(pvals.1)] <- 1
  signed.pvals.1[is.na(pvals.1)] <- 1
  
  rm(mat1,mat0,mat_all)
  
  #Rank p-values (ascending) for aggregation procedure
  struct_map$sign.pval <- rep(signed.pvals.1,times=as.vector(table(struct_map$L1)))
  struct_map <- struct_map %>% group_by(domainID,L1) %>% arrange(domainID,sign.pval) %>% ungroup() %>% data.table()
  #struct_map$rankL1 <- rleid(struct_map$L1)
  
  for(l in seq(L)){
    
    if(l==1){
      
      m.l <- length(pvals.1)
      
      pvals.l = pvals.1
      
      p.hat.l = est.p.hat(l=l,
                          D_prev = D_approx_prev,
                          FD_prev = FD_approx_prev,
                          pvals_l=pvals.l,
                          alpha=alpha)
      
      S.l <- which(pvals.l <= p.hat.l)
      
      if(teststat=="score"){
        LayerSumm[[l]]=list("score.l"=score.1,
                            "pvals.l"=pvals.l,
                            "S.l"=S.l,
                            "struct_map"=struct_map)
      } else{
        LayerSumm[[l]]=list("pvals.l"=pvals.l,
                            "S.l"=S.l,
                            "struct_map"=struct_map)
      }
      
    }else { #l>1
      
      #Update struct_map
      #First remove domains with fewer than 2^(l-1) L1 leaves
      struct_map <- struct_map[ , if (.N >= 2^(l-1) && uniqueN(L1) >= 2^(l-1)) .SD, by = .(domainID)]
      
      struct_map$rankL1 <- rleid(struct_map$L1)
      
      #Create layer l indices 
      Ll <- paste0("L",l)
      struct_map <-  struct_map[, (Ll) := {
        ## modify rle values
        x <- ceiling(rleid(rankL1) / 2^(l-1))
        n <- uniqueN(rankL1)
        if(n > 1 && between(n %% 2^(l-1),1,2^(l-1)-1)) {
          x[x == x[.N]] <- x[.N] - 1
        }
        x
      }, by = .(domainID)][, (Ll) := rleid(domainID, get(Ll))] #reassign groups by domain and Ll
      
      
      if(l==2){
        Llm1 <- paste0("rankL",l-1)
      } else{
        Llm1 <- paste0("L",l-1)
      }
      
      struct_map_tmp <- struct_map[, .(numSNP_lev = uniqueN(snpID),
                                       numL1_lev = uniqueN(rankL1),
                                       numLlm1_lev = uniqueN(get(Llm1))),
                                   by = .(domainID,get(Ll))]
      setnames(struct_map_tmp,old="get",new=Ll)
      
      #Number of layer 1 leaves used to create each layer l leaf
      numL1_leaves <- struct_map_tmp$numL1_lev
      
      numLlm1_leaves <- struct_map_tmp$numLlm1_lev
      
      #Subset layer 1 pvals to reflect current leaf indices
      #pvals.sub = pvals.lm1[unique(struct_map[,get(Llm1)])]
      
      #Get first obs. for each unique combo of (domainID, Ll, Llm1)
      pval_tmp <- struct_map %>%
        group_by(domainID,get(Ll),get(Llm1)) %>%
        slice(1)
      
      pvals.sub <- abs(pval_tmp$sign.pval)
      
      #Transform layer l-1 p-values before aggregating
      Q.sub = -2*log(pvals.sub)
      
      #Aggregate L(l-1) leaves according to partition defined on the lth layer
      Q.l <- as.vector(tapply(Q.sub,rep(seq_along(numLlm1_leaves),
                                        times=numLlm1_leaves),sum))
      
      m.l = length(Q.l)
      
      #Convert Q.l's to marginal p-values
      #Might be sub-uniform
      if(use.marg.p){
        pvals.l <- vapply(seq(m.l), function(x) pchisq(Q.l[x],df=2*numLlm1_leaves[x],
                                                       lower.tail = FALSE), numeric(1))
      } else{
        #Calculate conditional p-values depending on number of leaves merged from layer l-1
        pvals.l <- vapply(seq(m.l), function(x) ifelse(numLlm1_leaves[x]==2,
                                                       condChiSqProb2(Q.l[x],t_prev=-2*log(p.hats[l-1]),
                                                                      nu_prev_vec=c(2,2)),
                                                       condChiSqProb3(Q.l[x],t_prev=-2*log(p.hats[l-1]),
                                                                      nu_prev_vec=c(2,2,2))),
                          numeric(1))
      }
      
      
      #Obtain layer specific-threshold
      p.hat.l = est.p.hat(l=l,
                          D_prev = D_approx_prev,
                          FD_prev = FD_approx_prev,
                          pvals_l=pvals.l,
                          alpha=alpha)
      
      
      S.l <- which(pvals.l <= p.hat.l) #w.r.t. layer l indices
      
      
      LayerSumm[[l]]=list("Q.l"=Q.sub,
                          "pvals.l"=pvals.l,
                          "S.l"=S.l,
                          "struct_map"=struct_map)
      
    }
    
    
    #### Post-layer processing
    
    if(l==1){
      S = S.l
    } else{
      R.l <- unique(struct_map$L1[struct_map[,get(Ll)]%in%S.l])
      S = c(S,R.l)
    }
    
    #Save pvalues from previous layer
    #pvals.lm1 = pvals.l
    
    p.hats[l] = p.hat.l
    S.list[[l]] = S #discovered layer 1 leaves
    S.list.dc[[l]] = unique(struct_map_$domain[struct_map_$L1%in%S]) #vector of domain indices
    
    FD_approx_prev = FD_approx_prev + m.l*p.hats[l]
    D_approx_prev = length(S.l)
    
    #Update struct_map
    struct_map <- struct_map[!(struct_map$L1 %in% S), ]
    
  }
  
  
  return(list("params"=in.params,
              "p.hats"=p.hats,
              "S.list"=S.list,
              "S.list.dc"= S.list.dc,
              "LayerSumm"=LayerSumm))
}

### Simulation functions - each SNP is a leaf

SATET_sim_snp <- function(Gmat_case,Gmat_ctrl,
                          sim_map,
                          glm_input=NULL,
                          teststat="FET",
                          midp=TRUE,
                          score.cutoff=2,
                          use.marg.p=FALSE,
                          use.SPA.score=TRUE,
                          L=3,
                          alpha=0.05){
  
  require(purrr)
  require(data.table)
  require(dplyr)
  
  
  #Global objects
  p.hats <- rep(NA,L)
  S.l.list = vector("list",L)
  S.list=vector("list",L)
  S.list.dc=vector("list",L)
  LayerSumm = vector("list",L)
  
  FD_approx_prev = 0
  D_approx_prev = 0
  FD_true_prev = 0
  D_true_prev = 0
  
  #Simulation summary
  m <- ncol(Gmat_case)
  alt_snps <- sim_map$alt_snps
  alt_omni_snps <- sim_map$alt_omni
  null_snps <- setdiff(seq(m),c(alt_snps,alt_omni_snps))
  
  #Some processing
  if(class(Gmat_case)%in%c("CsparseMatrix","ngCMatrix")){
    mat1 <- unname(as(Gmat_case,"dgCMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgCMatrix"))
  } else{
    mat1 <- unname(as(Gmat_case,"dgTMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgTMatrix"))
  }
  
  mat_all <- rbind(mat1,mat0)
  
  N1 <- nrow(mat1)
  N0 <- nrow(mat0)
  N <- N1 + N0
  
  domain_snp_list <- sim_map$domain_snp_list
  n.snps <- sim_map$num_snps
  
  #Create DT that keeps track of leaf indices at each layer
  domain <- rep(seq_along(domain_snp_list), lengths(domain_snp_list))
  struct_map <- data.table(domainID = domain,
                           snpID = unlist(domain_snp_list),
                           key = "snpID")
  
  #Treat each snp as a bottom-layer leaf
  
  struct_map <- struct_map[, L1 := snpID]
  
  if(teststat=="FET"){
    #Get marginals - qualifying variants for cases and all
    
    case_colSums = Matrix::colSums(mat1)
    all_colSums = Matrix::colSums(mat_all)
    
    pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                    case_colSums = case_colSums,
                                    all_colSums = all_colSums,
                                    midp=midp)
    
  } else{
    #For score statistic
    leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
    
    if(!is.null(glm_input) && ncol(glm_input)>1){
      #With covariates
      if(use.SPA.score){
        score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                               pheno=glm_input[,1],cov=glm_input[,-1],
                                               minmac=1,Cutoff=score.cutoff)
      } else{
        score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                       pheno=glm_input[,1],cov=NULL,minmac=1)
      }
      
    } else{
      #Without covariates
      if(use.SPA.score){
        score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                               pheno=glm_input[,1],cov=NULL,
                                               minmac=1,Cutoff=score.cutoff)
      } else{
        score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                       pheno=glm_input[,1],cov=NULL,minmac=1)
      }
      
    }
    
    score.1 <- score.test$Tstat.sign
    pvals.1 <- score.test$p.value
  }
  
  
  #Copy object - this is not changed
  struct_map_ <- struct_map
  
  #If NaN, set to zero
  #(occurs when genotype vector is zero, so denominator of score stat is zero)
  pvals.1[is.na(pvals.1)] <- 1
  
  rm(mat1,mat0,mat_all)
  
  for(l in seq(L)){
    
    if(l==1){
      
      m.l <- length(pvals.1)
      
      pvals.l = pvals.1
      
      p.hat.l = est.p.hat(l=l,
                          D_prev = D_approx_prev,
                          FD_prev = FD_approx_prev,
                          pvals_l=pvals.l,
                          alpha=alpha)
      
      S.l <- which(pvals.l <= p.hat.l) 
      
      S.nul.lto1.len <- sum(S.l%in%null_snps)
      
      FD_true_l <- S.nul.lto1.len
      D_true_l <- length(S.l)
      
      FDP_true_l <- FD_true_l/max(D_true_l,1)
      
      if(teststat=="score"){
        LayerSumm[[1]]=list("pvals.l"=pvals.1,
                            "alt_snps"=alt_snps,
                            "alt_omni_snps"=alt_omni_snps,
                            "null_snps"=null_snps,
                            "FDP.l"=FDP_true_l,
                            "S.l"=S.l,
                            "struct_map"=struct_map_)
      } else{
        LayerSumm[[1]]=list("pvals.l"=pvals.1,
                            "alt_snps"=alt_snps,
                            "alt_omni_snps"=alt_omni_snps,
                            "null_snps"=null_snps,
                            "FDP.l"=FDP_true_l,
                            "S.l"=S.l,
                            "struct_map"=struct_map_)
      }
      
      
    }else { #l>1
      
      #Update leaf_indx_map
      #First remove domains with fewer than 2^(l-1) L1 leaves
      struct_map <- struct_map[ , if (.N >= 2^(l-1) && uniqueN(L1) >= 2^(l-1)) .SD, by = .(domainID)]
      
      Ll <- paste0("L",l)
      struct_map <-  struct_map[, (Ll) := {
        ## modify rle values
        x <- ceiling(rleid(L1) / 2^(l-1))
        n <- uniqueN(L1)
        if(n > 1 && between(n %% 2^(l-1),1,2^(l-1)-1)) {
          x[x == x[.N]] <- x[.N] - 1
        }
        x
      }, by = .(domainID)][, (Ll) := rleid(domainID, get(Ll))] #reassign groups by domain and Ll
      
      
      Llm1 <- paste0("L",l-1)
      struct_map_tmp <- struct_map[, .(numSNP_lev = uniqueN(snpID),
                                       numL1_lev = uniqueN(L1),
                                       numLlm1_lev = uniqueN(get(Llm1))),
                                   by = .(domainID,get(Ll))]
      setnames(struct_map_tmp,old="get",new=Ll)
      
      #Number of layer 1 leaves used to create each layer l leaf
      numL1_leaves <- struct_map_tmp$numL1_lev
      
      numLlm1_leaves <- struct_map_tmp$numLlm1_lev
      
      #Subset layer 1 pvals to reflect current leaf indices
      #Make sure to remove layer l-1 leaves from domains with less than 2^(l-1) L1 leaves
      pvals.sub = pvals.lm1[unique(struct_map[,get(Llm1)])]
      
      #Transform layer l-1 p-values before aggregating
      Q.sub = -2*log(pvals.sub)
      
      #Aggregate L(l-1) leaves according to partition defined on the lth layer
      Q.l <- as.vector(tapply(Q.sub,rep(seq_along(numLlm1_leaves),
                                        times=numLlm1_leaves),sum))
      
      m.l = length(Q.l)
      
      #Convert Q.l's to marginal p-values
      #Might be sub-uniform
      if(use.marg.p){
        pvals.l <- vapply(seq(m.l), function(x) pchisq(Q.l[x],df=2*numLlm1_leaves[x],
                                                       lower.tail = FALSE), numeric(1))
      } else{
        #Calculate conditional p-values depending on number of leaves merged from layer l-1
        pvals.l <- vapply(seq(m.l), function(x) ifelse(numLlm1_leaves[x]==2,
                                                       condChiSqProb2(Q.l[x],t_prev=-2*log(p.hats[l-1]),
                                                                      nu_prev_vec=c(2,2)),
                                                       condChiSqProb3(Q.l[x],t_prev=-2*log(p.hats[l-1]),
                                                                      nu_prev_vec=c(2,2,2))),
                          numeric(1))
      }
      
      
      #Obtain layer specific-threshold
      p.hat.l = est.p.hat(l=l,
                          D_prev = D_approx_prev,
                          FD_prev = FD_approx_prev,
                          pvals_l=pvals.l,
                          alpha=alpha)
      
      
      S.l <- which(pvals.l <= p.hat.l) #w.r.t. layer l indices
      
      S.nul.lto1.len <- sapply(S.l,
                               function(x) sum(unique(struct_map$L1[struct_map[,get(Ll)]==x])%in%null_snps))
      
      S.lto1.len <- sapply(S.l,
                           function(x) uniqueN(struct_map$L1[struct_map[,get(Ll)]==x]))
      
      if(length(S.nul.lto1.len)>0 & length(S.lto1.len)>0){
        FD_true_l <- FD_true_prev + sum(S.nul.lto1.len/pmax(S.lto1.len,1))
      } else{
        FD_true_l <- FD_true_prev
      }
      
      D_true_l <- D_true_prev + length(S.l)
      
      FDP_true_l <- FD_true_l/max(D_true_l,1)
      
      LayerSumm[[l]]=list("Q.l"=Q.sub,
                          "pvals.l"=pvals.l,
                          "FDP.l"=FDP_true_l, 
                          "S.l"=S.l,
                          "struct_map"=struct_map)
      
    }
    
    ##### Post-layer processing
    FD_true_prev <- FD_true_l
    D_true_prev <- D_true_l
    
    if(l==1){
      S = S.l
    } else{
      R.l <- unique(struct_map$L1[struct_map[,get(Ll)]%in%S.l])
      S = c(S,R.l)
    }
    
    #Save pvalues from previous layer
    pvals.lm1 = pvals.l
    
    p.hats[l] = p.hat.l
    S.list[[l]] = S #discovered layer 1 leaves
    S.list.dc[[l]] = unique(struct_map_$domain[struct_map_$L1%in%S]) #vector of domain indices
    
    FD_approx_prev = FD_approx_prev + m.l*p.hats[l]
    D_approx_prev = length(S.l)
    
    #Update struct_map
    struct_map <- struct_map[!(struct_map$L1 %in% S), ]
    
  }
  
  return(list("p.hats"=p.hats,
              "S.list"=S.list,
              "S.list.dc"=S.list.dc,
              "LayerSumm"=LayerSumm))
}

#### Simulation function - compare with domain collapsing

SATET_sim_dc <- function(Gmat_case,Gmat_ctrl,
                         sim_map,
                         glm_input=NULL,
                         teststat="FET",
                         midp=TRUE,
                         score.cutoff=2,
                         use.marg.p=TRUE,
                         use.SPA.score=TRUE,
                         alpha=0.05){
  
  require(purrr)
  require(data.table)
  require(dplyr)
  
  #Some processing
  if(class(Gmat_case)%in%c("CsparseMatrix","ngCMatrix")){
    mat1 <- unname(as(Gmat_case,"dgCMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgCMatrix"))
  } else{
    mat1 <- unname(as(Gmat_case,"dgTMatrix"))
    mat0 <- unname(as(Gmat_ctrl,"dgTMatrix"))
  }
  
  mat_all <- rbind(mat1,mat0)
  
  N1 <- nrow(mat1)
  N0 <- nrow(mat0)
  N <- N1 + N0
  
  domain_snp_list <- sim_map$domain_snp_list
  n.snps <- sim_map$num_snps
  
  #Create DT that keeps track of leaf indices at each layer
  domain <- rep(seq_along(domain_snp_list), lengths(domain_snp_list))
  struct_map <- data.table(domainID = domain,
                           snpID = unlist(domain_snp_list),
                           key = "snpID")
  
  struct_map <- struct_map[, L1 := domainID]
  
  leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
  
  if(teststat=="FET"){
    
    #leaf_mat_all <- create_leaf_attribute(mat_all,struct_map)
    leaf_mat_all@x <- ifelse(leaf_mat_all@x>0,1,0) #binarize
    
    case_colSums = Matrix::colSums(leaf_mat_all[1:N1,])
    all_colSums = Matrix::colSums(leaf_mat_all)
    
    pvals.1 <- calcFETpval_per_leaf(N1=N1,N0=N0,
                                    case_colSums = case_colSums,
                                    all_colSums = all_colSums,
                                    midp=midp)
    
  } else{
    
    
    if(!is.null(glm_input) && ncol(glm_input)>1){
      #With covariates
      if(use.SPA.score){
        score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                               pheno=glm_input[,1],cov=glm_input[,-1],
                                               minmac=1,Cutoff=score.cutoff)
      } else{
        score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                       pheno=glm_input[,1],cov=NULL,minmac=1)
      }
      
    } else{
      #Without covariates
      if(use.SPA.score){
        score.test <- ScoreTest_fastSPA_sparse(genomat=leaf_mat_all, #leaf_mat_all is an Nxm matrix of leaf attributes
                                               pheno=glm_input[,1],cov=NULL,
                                               minmac=1,Cutoff=score.cutoff)
      } else{
        score.test <- ScoreTest_sparse(genomat=leaf_mat_all,
                                       pheno=glm_input[,1],cov=NULL,minmac=1)
      }
      
    }
    
    score.1 <- score.test$Tstat.sign
    pvals.1 <- score.test$p.value
  }
  
  #If NaN, set to zero
  #(occurs when genotype vector is zero, so denominator of score stat is zero)
  pvals.1[is.na(pvals.1)] <- 1
  
  rm(mat1,mat0,mat_all)
  
  m.l <- length(pvals.1)
  
  p.hat.l = est.p.hat.dc(pvals_l=pvals.1,
                      alpha=alpha)
  
  S.1 <- which(pvals.1 <= p.hat.l) #Discoveries mapped to the leaf level
  
  
  #Simulation summary
  alt_snps <- sim_map$alt_snps
  alt_omni_snps <- sim_map$alt_omni
  
  alt_L1_leaves <- unique(struct_map$L1[struct_map$snpID%in%alt_snps])
  null_L1_leaves <- setdiff(seq(m.l),alt_L1_leaves)
  
  if(teststat=="score"){
    LayerSumm =list("score.l"=score.1,
                    "pvals.l"=pvals.1,
                    "alt_snps"=alt_snps,
                    "alt_omni_snps"=alt_omni_snps,
                    "null_snps"=setdiff(seq(n.snps),c(alt_snps,alt_omni_snps)),
                    "alt_leaves"=alt_L1_leaves,
                    "null_leaves"=null_L1_leaves,
                    "S.l"=S.1,
                    "struct_map"=struct_map)
  } else{
    LayerSumm =list("pvals.l"=pvals.1,
                    "alt_snps"=alt_snps,
                    "alt_omni_snps"=alt_omni_snps,
                    "null_snps"=setdiff(seq(n.snps),c(alt_snps,alt_omni_snps)),
                    "alt_leaves"=alt_L1_leaves,
                    "null_leaves"=null_L1_leaves,
                    "S.l"=S.1,
                    "struct_map"=struct_map)
  }
  
  
  return(list("p.hats"=p.hat.l,
              "S.list"=S.1,
              "S.list.dc"=S.1,
              "LayerSumm"=LayerSumm))
}




############# AUXILIARY FUNCTIONS #################
splitNoOverlap <- function(x,n){
  numOfVectors <- length(x)%/%n
  elementsPerVector <- c(rep(n,numOfVectors-1),n+length(x) %% n)
  elemDistPerVector <- rep(1:numOfVectors,elementsPerVector)
  split(x,factor(elemDistPerVector))
}

condChiSqProb2 <- function(t,t_prev,nu_prev_vec){
  require(pracma)
  if(length(nu_prev_vec)!=2){
    return(0)
  } else{
    if(t<=0){
      prob = 1
    } else if(t>=0 & t<=t_prev){
      integrand <- function(y) pchisq(t-y,df=nu_prev_vec[1])*dchisq(y,df=nu_prev_vec[2])
      I <- integral(integrand, 0, t)
      prob = 1-I/(pchisq(t_prev,df=nu_prev_vec[1])*pchisq(t_prev,df=nu_prev_vec[2]))
    } else if(t>t_prev & t<=2*t_prev){
      integrand <- function(y) (pchisq(t_prev,df=nu_prev_vec[1])-
                                  pchisq(t-y,df=nu_prev_vec[1]))*dchisq(y,df=nu_prev_vec[2])
      I <- integral(integrand, t-t_prev,t_prev)
      prob = I/(pchisq(t_prev,df=nu_prev_vec[1])*pchisq(t_prev,df=nu_prev_vec[2]))
    } else{
      prob = 0
    }
    return(prob)
  }
}

condChiSqProb3 <- function(t,t_prev,nu_prev_vec){
  require(pracma)
  if(length(nu_prev_vec)!=3){
    return(0)
  } else{
    if(t<=0){
      prob = 1
    } else if(t>0 & t<=t_prev){
      integrand <- function(x,y) pchisq(t-(x+y),df=nu_prev_vec[1])*dchisq(x,df=nu_prev_vec[2])*dchisq(y,df=nu_prev_vec[3])
      xmin <- 0; xmax <- t
      ymin <- 0; ymax <- function(x) t-x
      I <- integral2(integrand, xmin, xmax, ymin, ymax)
      prob = 1-I$Q/(pchisq(t_prev,df=nu_prev_vec[1])*pchisq(t_prev,df=nu_prev_vec[2])*pchisq(t_prev,df=nu_prev_vec[3]))
    } else if(t>t_prev & t<=2*t_prev){
      #Need to integrate over two 3D regions
      integrand <- function(x,y) (pchisq(t_prev,df=nu_prev_vec[1])-
                                    pchisq(t-(x+y),df=nu_prev_vec[1]))*dchisq(x,df=nu_prev_vec[2])*dchisq(y,df=nu_prev_vec[3])
      xmin1 <- 0; xmax1 <- t-t_prev
      ymin1 <- function(x) t-t_prev-x; ymax1 <- t_prev
      I1 <- integral2(integrand, xmin1, xmax1, ymin1, ymax1)
      
      xmin2 <- t-t_prev; xmax2 <- t_prev
      ymin2 <- 0; ymax2 <- function(x) t-x
      I2 <- integral2(integrand, xmin2, xmax2, ymin2, ymax2)
      
      prob = (I1$Q+I2$Q)/(pchisq(t_prev,df=nu_prev_vec[1])*pchisq(t_prev,df=nu_prev_vec[2])*pchisq(t_prev,df=nu_prev_vec[3]))
    } else if(t>2*t_prev & t<=3*t_prev){
      integrand <- function(x,y) (pchisq(t_prev,df=nu_prev_vec[1])-
                                    pchisq(t-(x+y),df=nu_prev_vec[1]))*dchisq(x,df=nu_prev_vec[2])*dchisq(y,df=nu_prev_vec[3])
      xmin <- t-2*t_prev; xmax <- t_prev
      ymin <- function(x) t-x; ymax <- t_prev
      I <- integral2(integrand, xmin, xmax, ymin, ymax)
      prob = I$Q/(pchisq(t_prev,df=nu_prev_vec[1])*pchisq(t_prev,df=nu_prev_vec[2])*pchisq(t_prev,df=nu_prev_vec[3]))
    } else{
      prob = 0
    }
    return(prob)
  }
}


est.p.hat <- function(l,D_prev,FD_prev,pvals_l,alpha){
  
  ##print(paste("l:",l))
  
  m.l = length(pvals_l)
  
  #Threshold very small p-values based on constant
  p.m = 1/(m.l*sqrt(log(m.l)))
  
  filter <- which(pvals_l<=p.m)
  
  if(length(filter)>0){
    p.vec = sort(pvals_l[-filter],decreasing = FALSE)
  }else {
    p.vec = sort(pvals_l,decreasing = FALSE)
  }
  
  p.indx = 0
  emp.fdr = 0
  
  while(all(emp.fdr <= alpha, p.indx < length(p.vec))){
    p.indx = p.indx + 1
      
    fdr.num = FD_prev + m.l*p.vec[p.indx]
      
    S.l <- which(pvals_l <= p.vec[p.indx]) #w.r.t. layer l indices
    fdr.denom = D_prev + length(S.l)
    
    emp.fdr = fdr.num/max(fdr.denom,1)
    
    #print(paste("p.indx:",p.indx))
    #print(paste("p.val:",p.vec[p.indx]))
    #print(paste("fdr.num:",fdr.num))
    #print(paste("fdr.denom:",fdr.denom))
    #print(paste("emp.fdr:",emp.fdr))
  }
  
  p.hat = ifelse(is.na(p.vec[p.indx-1]) || length(p.vec[p.indx-1])==0,
                 p.m,
                 p.vec[p.indx-1])
  
  return(p.hat)
}

est.p.hat.dc <- function(pvals_l,alpha){
  
  ##print(paste("l:",l))
  
  m.l = length(pvals_l)
  
  #Threshold very small p-values based on constant
  p.m = 1/(m.l*sqrt(log(m.l)))
  
  filter <- which(pvals_l<=p.m)

  if(length(filter)>0){
    p.vec = sort(pvals_l[-filter],decreasing = FALSE)
  }else {
    p.vec = sort(pvals_l,decreasing = FALSE)
  }
  
  p.indx = 0
  emp.fdr = 0
  
  while(all(emp.fdr <= alpha, p.indx < length(p.vec))){
    p.indx = p.indx + 1
    
    fdr.num = m.l*p.vec[p.indx]
    
    fdr.denom = sum(pvals_l <= p.vec[p.indx])
    emp.fdr = fdr.num/max(fdr.denom,1)
    
    #print(paste("p.indx:",p.indx))
    #print(paste("p.val:",p.vec[p.indx]))
    #print(paste("fdr.num:",fdr.num))
    #print(paste("fdr.denom:",fdr.denom))
    #print(paste("emp.fdr:",emp.fdr))
  }
  
  p.hat = ifelse(is.na(p.vec[p.indx-1]) || length(p.vec[p.indx-1])==0,
                 p.m,
                 p.vec[p.indx-1])
  
  return(p.hat)}


####### Functions to compute test-statistics

fisher.exact.test <- function(z,midp=TRUE){
  
  x <- z[1]
  sampTot <- z[2]
  pop1Tot <- z[3]
  pop2Tot <- z[4]
  
  lo <- max(0L, sampTot - pop2Tot)
  hi <- min(sampTot, pop1Tot)
  
  support <- lo : hi
  out <- dhyper(support, pop1Tot, pop2Tot, sampTot)
  
  if(midp){
    #mid p-val with minimum likelihood method
    return(sum(out[out < out[x - lo + 1]]) + sum(out[out==out[x-lo+1]])/2)
  } else{
    #minimum likelihood method
    return(sum(out[out <= out[x - lo + 1]]))
  }
}

calcFETpval_per_leaf <- function(N1,N0,case_colSums,all_colSums,
                                 midp=TRUE){
  
  #Order of input is: sample hit, sample size, pop1 size, pop2 size
  cont_tab_summ2x2 <- unname(cbind(case_colSums,all_colSums,N1,N0))
  
  # Apply row-wise
  FET_pvals <- apply(cont_tab_summ2x2, 1, 
                     function(z) fisher.exact.test(z,midp = midp))
  
  #Vector of length m.l, where m.l is number of leaves/hypotheses at layer l
  return(FET_pvals)
  
}

calc_pval_per_leaf2 <- function(case_Q_ind,ctrl_Q_ind,colSums_Q){
  
  #row margins
  row_margins <- c(unlist(lapply(case_Q_ind,length)),
                   unlist(lapply(ctrl_Q_ind,length)))
  
  cont_arr_summ <- array(NA,c(length(row_margins),2,ncol(colSums_Q)))
  
  for(i in seq(ncol(colSums_Q))){
    cont_arr_summ[,,i] <- cbind(colSums_Q[,i],row_margins-colSums_Q[,i])
  }
  
  chisq_pvals <- apply(cont_arr_summ, 3,
                     function(z) chisq.test(z)$p.value)
  # FET2_pvals <- apply(cont_arr_summ, 3,
  #                    function(z) fisher.test(z)$p.value)
  
  #Vector of length m.l, where m.l is number of leaves/hypotheses at layer l
  return(chisq_pvals)
  #return(FET2_pvals)
  
}



# CreateSNPLeafSet <- function(x,cutoff){
#   
#   groupVec = rep(NA,length(x))
#   group = 1
#   runSum = 0
#   for(i in seq_along(x)){
#     runSum <- runSum + x[i]
#     groupVec[i] <- group
#     if(runSum >= cutoff){
#       group <- group + 1
#       runSum = 0
#     }
#   }
#   
#   return(groupVec)
# }

sum_reset_at = function(val_col, threshold, include.equals = TRUE) {
  if (include.equals) {
    purrr::accumulate({{val_col}}, ~if_else(.x>=threshold , .y, .x+.y))
  } else {
    purrr::accumulate({{val_col}}, ~if_else(.x>threshold , .y, .x+.y))
  }
}

compute_counts <- function(thresh, ID, domain_ends) {
  
  seen_ids <- NULL
  count <- 0L
  adjust_count <- function(id, domain_end) {
    if (!(id %in% seen_ids)) {
      seen_ids <<- c(seen_ids,id)
      count <<- count + 1L
    }
    
    if (domain_end | uniqueN(seen_ids) >= thresh) {
      count <- count # copy enclosed value locally
      seen_ids <<- NULL
      count <<- 0L
    }
    
    count
  }
  
  unlist(Map(adjust_count, ID, domain_ends))
}


create_leaf_attribute <- function(snp_mat,snp_leaf_map){
  
  snp_mat@Dimnames <- list(NULL,NULL)
  
  #Convert mat to dgTMatrix if not already
  if(class(snp_mat)!="dgTMatrix"){
    snp_mat <- as(snp_mat, "dgTMatrix") 
  }
  
  snp_mat2 <- snp_mat #copy object
  
  #Replace column indices with new set of indices 
  #Make sure initial indices start with zero
  snp_mat2@j <- as.integer(snp_leaf_map$L1[snp_mat@j+1]-1)
  #Correct dimensions of new matrix
  
  snp_mat2@Dim <- as.integer(c(nrow(snp_mat2),length(unique(snp_leaf_map$L1))))
  
  #Convert to dgCMatrix
  y <- as(snp_mat2,"dgCMatrix")
  return(y)
}


############# SIMULATION FUNCTIONS ############

run_DGP1_sing <- function(N1,
                          sampsize.factors,
                          sim_map,
                          params,
                          midp=TRUE,
                          use.marg.p=FALSE,
                          covars=FALSE,
                          alpha,L,seed){

  set.seed(seed)
  

  Sim.mat = array(NA,c(1,length(sampsize.factors),9,L,7)) #1 iteration, prevalence, 3 teststats x 2 collapsing strategies, L layers, and 6 performance metrics

  
  ### Generate Sample Data ###
  p0 = 1/100
  
  for(i in seq_along(sampsize.factors)){
    
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP1(Npop=ceiling(3*N1/p0),
                    N0 = N0,N1 = N1,
                    sim_map = sim_map,
                    covars=covars,
                    params = params,
                    p0 = p0)
    
    
      res_FET_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                         Gmat_ctrl=attrDat$geno_ctrl_matrix,
                         sim_map=sim_map,
                         glm_input = attrDat$nulldata,
                         teststat = "FET",
                         use.marg.p = use.marg.p,
                         alpha=alpha)
      
      #Run algorithm with method DD and score statistic (logistic regression)
      res_score_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                           Gmat_ctrl=attrDat$geno_ctrl_matrix,
                           sim_map=sim_map,
                           glm_input = attrDat$nulldata,
                           teststat = "score",
                           score.cutoff = params$score.cutoff,
                           use.SPA.score = FALSE,
                           use.marg.p = use.marg.p,
                           alpha=alpha)
      
      res_score_wSPA_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                                Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                sim_map=sim_map,
                                glm_input = attrDat$nulldata,
                                teststat = "score",
                                score.cutoff = params$score.cutoff,
                                use.SPA.score = TRUE,
                                use.marg.p = use.marg.p,
                                alpha=alpha)
      
      #Run algorithm with method DD and FET statistic
      res_FET =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                         Gmat_ctrl=attrDat$geno_ctrl_matrix,
                         sim_map=sim_map,
                         glm_input = attrDat$nulldata,
                         teststat = "FET",
                         use.marg.p = use.marg.p,
                         L=L,
                         alpha=alpha)
      
      #Run algorithm with method DD and score statistic (logistic regression)
      res_score =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                           Gmat_ctrl=attrDat$geno_ctrl_matrix,
                           sim_map=sim_map,
                           glm_input = attrDat$nulldata,
                           teststat = "score",
                           score.cutoff = params$score.cutoff,
                           use.SPA.score = FALSE,
                           use.marg.p = use.marg.p,
                           L=L,
                           alpha=alpha)
      
      res_score_wSPA =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                                 Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                 sim_map=sim_map,
                                 glm_input = attrDat$nulldata,
                                 teststat = "score",
                                 score.cutoff = params$score.cutoff,
                                 use.SPA.score = TRUE,
                                use.marg.p = use.marg.p,
                                 L=L,
                                 alpha=alpha)
    
    #Performance
    
    #Ground truth (leaves)
    m = sim_map$num_snps
    nonnulls = sim_map$alt_snps
    nonnulls_omni = sim_map$alt_omni_snps
    nulls = sim_map$null_snps
    
    #Ground truth - domain collapsing
    m_dc = length(sim_map$domain_snp_list)
    nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
    nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
    
    for (l in 1:L){
      
        #Method FET + DD
        Sim.mat[1,i,1,l,1] = sum(res_FET$S.list[[l]]%in%nonnulls) #TP
        Sim.mat[1,i,1,l,2] = sum(res_FET$S.list[[l]]%in%nonnulls_omni) #TP omni
        Sim.mat[1,i,1,l,3] = sum(res_FET$S.list[[l]]%in%nulls) #FP
        Sim.mat[1,i,1,l,4] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nulls) #TN
        Sim.mat[1,i,1,l,5] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls) #FN
        Sim.mat[1,i,1,l,6] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls_omni) #FN omni
        Sim.mat[1,i,1,l,7] = res_FET$LayerSumm[[l]]$FDP.l
        
        #Method score + noSPA + DD
        Sim.mat[1,i,2,l,1] = sum(res_score$S.list[[l]]%in%nonnulls) #TP
        Sim.mat[1,i,2,l,2] = sum(res_score$S.list[[l]]%in%nonnulls_omni) #TP omni
        Sim.mat[1,i,2,l,3] = sum(res_score$S.list[[l]]%in%nulls) #FP
        Sim.mat[1,i,2,l,4] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nulls) #TN
        Sim.mat[1,i,2,l,5] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls) #FN
        Sim.mat[1,i,2,l,6] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls_omni) #FN omni
        Sim.mat[1,i,2,l,7] = res_score$LayerSumm[[l]]$FDP.l
        
        #Method scoreSPA + DD
        Sim.mat[1,i,3,l,1] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls) #TP
        Sim.mat[1,i,3,l,2] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls_omni) #TP omni
        Sim.mat[1,i,3,l,3] = sum(res_score_wSPA$S.list[[l]]%in%nulls) #FP
        Sim.mat[1,i,3,l,4] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nulls) #TN
        Sim.mat[1,i,3,l,5] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls) #FN
        Sim.mat[1,i,3,l,6] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls_omni) #FN omni
        Sim.mat[1,i,3,l,7] = res_score_wSPA$LayerSumm[[l]]$FDP.l
        
        #Method FET + DD + domain collapsing (Power comparison only)
        Sim.mat[1,i,4,l,1] = sum(res_FET$S.list.dc[[l]]%in%nonnulls_dc) #TP
        Sim.mat[1,i,4,l,2] = sum(res_FET$S.list.dc[[l]]%in%nulls_dc) #FP
        Sim.mat[1,i,4,l,3] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nulls_dc) #TN
        Sim.mat[1,i,4,l,4] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nonnulls_dc) #FN
        
        #Method score + DD + domain collapsing (Power comparison only)
        Sim.mat[1,i,5,l,1] = sum(res_score$S.list.dc[[l]]%in%nonnulls_dc) #TP
        Sim.mat[1,i,5,l,2] = sum(res_score$S.list.dc[[l]]%in%nulls_dc) #FP
        Sim.mat[1,i,5,l,3] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nulls_dc) #TN
        Sim.mat[1,i,5,l,4] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nonnulls_dc) #FN
        
        #Method score + DD + domain collapsing (Power comparison only)
        Sim.mat[1,i,6,l,1] = sum(res_score_wSPA$S.list.dc[[l]]%in%nonnulls_dc) #TP
        Sim.mat[1,i,6,l,2] = sum(res_score_wSPA$S.list.dc[[l]]%in%nulls_dc) #FP
        Sim.mat[1,i,6,l,3] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nulls_dc) #TN
        Sim.mat[1,i,6,l,4] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nonnulls_dc) #FN
    }
    
    #Method FET + domain collapsing
    Sim.mat[1,i,7,1,1] = sum(res_FET_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,7,1,2] = sum(res_FET_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,7,1,3] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,7,1,4] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nonnulls_dc) #FN
 
    #Method score + domain collapsing
    Sim.mat[1,i,8,1,1] = sum(res_score_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,8,1,2] = sum(res_score_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,8,1,3] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,8,1,4] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nonnulls_dc) #FN

    #Method score + domain collapsing
    Sim.mat[1,i,9,1,1] = sum(res_score_wSPA_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,9,1,2] = sum(res_score_wSPA_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,9,1,3] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,9,1,4] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nonnulls_dc) #FN

  }
  
  return(Sim.mat)
}

run_DGP1_sing_conf <- function(N1,
                          sampsize.factors,
                          sim_map,
                          params,
                          midp=TRUE,
                          use.marg.p=FALSE,
                          covars=FALSE,
                          alpha,L,seed){
  
  set.seed(seed)
  
  
  Sim.mat = array(NA,c(1,length(sampsize.factors),9,L,7)) #1 iteration, prevalence, 3 teststats x 2 collapsing strategies, L layers, and 6 performance metrics
  
  
  ### Generate Sample Data ###
  p0 = 1/100
  
  for(i in seq_along(sampsize.factors)){
    
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP1(Npop=ceiling(3*N1/p0),
                    N0 = N0,N1 = N1,
                    sim_map = sim_map,
                    covars=covars,
                    params = params,
                    p0 = p0)
    
    
    res_FET_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                             Gmat_ctrl=attrDat$geno_ctrl_matrix,
                             sim_map=sim_map,
                             glm_input = attrDat$nulldata,
                             teststat = "FET",
                             use.marg.p = use.marg.p,
                             alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                               Gmat_ctrl=attrDat$geno_ctrl_matrix,
                               sim_map=sim_map,
                               glm_input = attrDat$nulldata,
                               teststat = "score",
                               score.cutoff = params$score.cutoff,
                               use.SPA.score = FALSE,
                               use.marg.p = use.marg.p,
                               alpha=alpha)
    
    res_score_wSPA_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                                    Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                    sim_map=sim_map,
                                    glm_input = attrDat$nulldata,
                                    teststat = "score",
                                    score.cutoff = params$score.cutoff,
                                    use.SPA.score = TRUE,
                                    use.marg.p = use.marg.p,
                                    alpha=alpha)
    
    #Run algorithm with method DD and FET statistic
    res_FET =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                           Gmat_ctrl=attrDat$geno_ctrl_matrix,
                           sim_map=sim_map,
                           glm_input = attrDat$nulldata,
                           teststat = "FET",
                           use.marg.p = use.marg.p,
                           L=L,
                           alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                             Gmat_ctrl=attrDat$geno_ctrl_matrix,
                             sim_map=sim_map,
                             glm_input = attrDat$nulldata,
                             teststat = "score",
                             score.cutoff = params$score.cutoff,
                             use.SPA.score = FALSE,
                             use.marg.p = use.marg.p,
                             L=L,
                             alpha=alpha)
    
    res_score_wSPA =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                                  Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                  sim_map=sim_map,
                                  glm_input = attrDat$nulldata,
                                  teststat = "score",
                                  score.cutoff = params$score.cutoff,
                                  use.SPA.score = TRUE,
                                  use.marg.p = use.marg.p,
                                  L=L,
                                  alpha=alpha)
    
    #Performance
    
    #Ground truth (leaves)
    m = sim_map$num_snps
    nonnulls = sim_map$alt_snps
    nonnulls_omni = sim_map$alt_omni_snps
    nulls = sim_map$null_snps
    
    #Ground truth - domain collapsing
    m_dc = length(sim_map$domain_snp_list)
    nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
    nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
    
    for (l in 1:L){
      
      #Method FET + DD
      Sim.mat[1,i,1,l,1] = sum(res_FET$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,1,l,2] = sum(res_FET$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,1,l,3] = sum(res_FET$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,1,l,4] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,1,l,5] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,1,l,6] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls_omni) #FN omni
      Sim.mat[1,i,1,l,7] = res_FET$LayerSumm[[l]]$FDP.l
      
      #Method score + noSPA + DD
      Sim.mat[1,i,2,l,1] = sum(res_score$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,2,l,2] = sum(res_score$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,2,l,3] = sum(res_score$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,2,l,4] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,2,l,5] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,2,l,6] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls_omni) #FN omni
      Sim.mat[1,i,2,l,7] = res_score$LayerSumm[[l]]$FDP.l
      
      #Method scoreSPA + DD
      Sim.mat[1,i,3,l,1] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,3,l,2] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,3,l,3] = sum(res_score_wSPA$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,3,l,4] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,3,l,5] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,3,l,6] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls_omni) #FN omni
      Sim.mat[1,i,3,l,7] = res_score_wSPA$LayerSumm[[l]]$FDP.l
      
      #Method FET + DD + domain collapsing (Power comparison only)
      Sim.mat[1,i,4,l,1] = sum(res_FET$S.list.dc[[l]]%in%nonnulls_dc) #TP
      Sim.mat[1,i,4,l,2] = sum(res_FET$S.list.dc[[l]]%in%nulls_dc) #FP
      Sim.mat[1,i,4,l,3] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nulls_dc) #TN
      Sim.mat[1,i,4,l,4] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nonnulls_dc) #FN
      
      #Method score + DD + domain collapsing (Power comparison only)
      Sim.mat[1,i,5,l,1] = sum(res_score$S.list.dc[[l]]%in%nonnulls_dc) #TP
      Sim.mat[1,i,5,l,2] = sum(res_score$S.list.dc[[l]]%in%nulls_dc) #FP
      Sim.mat[1,i,5,l,3] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nulls_dc) #TN
      Sim.mat[1,i,5,l,4] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nonnulls_dc) #FN
      
      #Method score + DD + domain collapsing (Power comparison only)
      Sim.mat[1,i,6,l,1] = sum(res_score_wSPA$S.list.dc[[l]]%in%nonnulls_dc) #TP
      Sim.mat[1,i,6,l,2] = sum(res_score_wSPA$S.list.dc[[l]]%in%nulls_dc) #FP
      Sim.mat[1,i,6,l,3] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nulls_dc) #TN
      Sim.mat[1,i,6,l,4] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nonnulls_dc) #FN
    }
    
    #Method FET + domain collapsing
    Sim.mat[1,i,7,1,1] = sum(res_FET_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,7,1,2] = sum(res_FET_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,7,1,3] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,7,1,4] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nonnulls_dc) #FN
    
    #Method score + domain collapsing
    Sim.mat[1,i,8,1,1] = sum(res_score_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,8,1,2] = sum(res_score_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,8,1,3] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,8,1,4] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nonnulls_dc) #FN
    
    #Method score + domain collapsing
    Sim.mat[1,i,9,1,1] = sum(res_score_wSPA_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,9,1,2] = sum(res_score_wSPA_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,9,1,3] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,9,1,4] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nonnulls_dc) #FN
    
  }
  
  return(Sim.mat)
}


run_DGP2_sing <- function(N1,
                          sampsize.factors,
                          sim_map,
                          params,
                          midp=TRUE,
                          use.marg.p=FALSE,
                          covars=FALSE,
                          alpha,L,seed){
  
  set.seed(seed)
  
  Sim.mat = array(NA,c(1,length(sampsize.factors),9,L,7)) #1 iteration, prevalence, 3 teststats x 2 collapsing strategies, L layers, and 6 performance metrics

  
  ### Generate Sample Data ###
  p0 = 1/100
  
  for(i in seq_along(sampsize.factors)){
    
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP2(Npop=ceiling(3*N1/p0),
                    N0 = N0,N1 = N1,
                    sim_map = sim_map,
                    covars=covars,
                    params = params,
                    p0 = p0)
    
    res_FET_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                             Gmat_ctrl=attrDat$geno_ctrl_matrix,
                             sim_map=sim_map,
                             glm_input = attrDat$nulldata,
                             teststat = "FET",
                             use.marg.p = use.marg.p,
                             alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                               Gmat_ctrl=attrDat$geno_ctrl_matrix,
                               sim_map=sim_map,
                               glm_input = attrDat$nulldata,
                               teststat = "score",
                               score.cutoff = params$score.cutoff,
                               use.SPA.score = FALSE,
                               use.marg.p = use.marg.p,
                               alpha=alpha)
    
    res_score_wSPA_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                                    Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                    sim_map=sim_map,
                                    glm_input = attrDat$nulldata,
                                    teststat = "score",
                                    score.cutoff = params$score.cutoff,
                                    use.SPA.score = TRUE,
                                    use.marg.p = use.marg.p,
                                    alpha=alpha)
    
    
    #Run algorithm with method DD and FET statistic
    res_FET =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                       Gmat_ctrl=attrDat$geno_ctrl_matrix,
                       sim_map=sim_map,
                       glm_input = attrDat$nulldata,
                       teststat = "FET",
                       use.marg.p = use.marg.p,
                       L=L,
                       alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                         Gmat_ctrl=attrDat$geno_ctrl_matrix,
                         sim_map=sim_map,
                         glm_input = attrDat$nulldata,
                         teststat = "score",
                         score.cutoff = params$score.cutoff,
                         use.SPA.score = FALSE,
                         use.marg.p = use.marg.p,
                         L=L,
                         alpha=alpha)
    
    res_score_wSPA =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                              Gmat_ctrl=attrDat$geno_ctrl_matrix,
                              sim_map=sim_map,
                              glm_input = attrDat$nulldata,
                              teststat = "score",
                              score.cutoff = params$score.cutoff,
                              use.SPA.score = TRUE,
                              use.marg.p = use.marg.p,
                              L=L,
                              alpha=alpha)
    
    #Performance
    
    #Ground truth (leaves)
    m = sim_map$num_snps
    nonnulls = sim_map$alt_snps
    nonnulls_omni = sim_map$alt_omni_snps
    nulls = sim_map$null_snps
    
    #Ground truth - domain collapsing
    m_dc = length(sim_map$domain_snp_list)
    nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
    nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
    
    for (l in 1:L){
        
        #Method FET + DD
        Sim.mat[1,i,1,l,1] = sum(res_FET$S.list[[l]]%in%nonnulls) #TP
        Sim.mat[1,i,1,l,2] = sum(res_FET$S.list[[l]]%in%nonnulls_omni) #TP omni
        Sim.mat[1,i,1,l,3] = sum(res_FET$S.list[[l]]%in%nulls) #FP
        Sim.mat[1,i,1,l,4] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nulls) #TN
        Sim.mat[1,i,1,l,5] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls) #FN
        Sim.mat[1,i,1,l,6] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls_omni) #FN omni
        Sim.mat[1,i,1,l,7] = res_FET$LayerSumm[[l]]$FDP.l
        
        #Method score + noSPA + DD
        Sim.mat[1,i,2,l,1] = sum(res_score$S.list[[l]]%in%nonnulls) #TP
        Sim.mat[1,i,2,l,2] = sum(res_score$S.list[[l]]%in%nonnulls_omni) #TP omni
        Sim.mat[1,i,2,l,3] = sum(res_score$S.list[[l]]%in%nulls) #FP
        Sim.mat[1,i,2,l,4] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nulls) #TN
        Sim.mat[1,i,2,l,5] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls) #FN
        Sim.mat[1,i,2,l,6] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls_omni) #FN omni
        Sim.mat[1,i,2,l,7] = res_score$LayerSumm[[l]]$FDP.l
        
        #Method scoreSPA + DD
        Sim.mat[1,i,3,l,1] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls) #TP
        Sim.mat[1,i,3,l,2] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls_omni) #TP omni
        Sim.mat[1,i,3,l,3] = sum(res_score_wSPA$S.list[[l]]%in%nulls) #FP
        Sim.mat[1,i,3,l,4] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nulls) #TN
        Sim.mat[1,i,3,l,5] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls) #FN
        Sim.mat[1,i,3,l,6] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls_omni) #FN omni
        Sim.mat[1,i,3,l,7] = res_score_wSPA$LayerSumm[[l]]$FDP.l
        
        #Method FET + DD + domain collapsing (Power comparison only)
        Sim.mat[1,i,4,l,1] = sum(res_FET$S.list.dc[[l]]%in%nonnulls_dc) #TP
        Sim.mat[1,i,4,l,2] = sum(res_FET$S.list.dc[[l]]%in%nulls_dc) #FP
        Sim.mat[1,i,4,l,3] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nulls_dc) #TN
        Sim.mat[1,i,4,l,4] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nonnulls_dc) #FN
        
        #Method score + DD + domain collapsing (Power comparison only)
        Sim.mat[1,i,5,l,1] = sum(res_score$S.list.dc[[l]]%in%nonnulls_dc) #TP
        Sim.mat[1,i,5,l,2] = sum(res_score$S.list.dc[[l]]%in%nulls_dc) #FP
        Sim.mat[1,i,5,l,3] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nulls_dc) #TN
        Sim.mat[1,i,5,l,4] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nonnulls_dc) #FN
        
        #Method score + DD + domain collapsing (Power comparison only)
        Sim.mat[1,i,6,l,1] = sum(res_score_wSPA$S.list.dc[[l]]%in%nonnulls_dc) #TP
        Sim.mat[1,i,6,l,2] = sum(res_score_wSPA$S.list.dc[[l]]%in%nulls_dc) #FP
        Sim.mat[1,i,6,l,3] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nulls_dc) #TN
        Sim.mat[1,i,6,l,4] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nonnulls_dc) #FN
      }
      
      #Method FET + domain collapsing
      Sim.mat[1,i,7,1,1] = sum(res_FET_dc$S.list.dc%in%nonnulls_dc) #TP
      Sim.mat[1,i,7,1,2] = sum(res_FET_dc$S.list.dc%in%nulls_dc) #FP
      Sim.mat[1,i,7,1,3] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nulls_dc) #TN
      Sim.mat[1,i,7,1,4] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nonnulls_dc) #FN
      
      #Method score + domain collapsing
      Sim.mat[1,i,8,1,1] = sum(res_score_dc$S.list.dc%in%nonnulls_dc) #TP
      Sim.mat[1,i,8,1,2] = sum(res_score_dc$S.list.dc%in%nulls_dc) #FP
      Sim.mat[1,i,8,1,3] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nulls_dc) #TN
      Sim.mat[1,i,8,1,4] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nonnulls_dc) #FN
      
      #Method score + domain collapsing
      Sim.mat[1,i,9,1,1] = sum(res_score_wSPA_dc$S.list.dc%in%nonnulls_dc) #TP
      Sim.mat[1,i,9,1,2] = sum(res_score_wSPA_dc$S.list.dc%in%nulls_dc) #FP
      Sim.mat[1,i,9,1,3] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nulls_dc) #TN
      Sim.mat[1,i,9,1,4] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nonnulls_dc) #FN
    
  }
  
  return(Sim.mat)
}

run_DGP2_sing_conf <- function(N1,
                          sampsize.factors,
                          sim_map,
                          params,
                          midp=TRUE,
                          use.marg.p=FALSE,
                          covars=FALSE,
                          alpha,L,seed){
  
  set.seed(seed)
  
  Sim.mat = array(NA,c(1,length(sampsize.factors),9,L,7)) #1 iteration, prevalence, 3 teststats x 2 collapsing strategies, L layers, and 6 performance metrics
  
  
  ### Generate Sample Data ###
  p0 = 1/100
  
  for(i in seq_along(sampsize.factors)){
    
    ## Generate case-control samples
    N0 <- ceiling(N1*sampsize.factors[i])
    
    ## Create data
    attrDat <- DGP2(Npop=ceiling(3*N1/p0),
                    N0 = N0,N1 = N1,
                    sim_map = sim_map,
                    covars=covars,
                    params = params,
                    p0 = p0)
    
    res_FET_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                             Gmat_ctrl=attrDat$geno_ctrl_matrix,
                             sim_map=sim_map,
                             glm_input = attrDat$nulldata,
                             teststat = "FET",
                             use.marg.p = use.marg.p,
                             alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                               Gmat_ctrl=attrDat$geno_ctrl_matrix,
                               sim_map=sim_map,
                               glm_input = attrDat$nulldata,
                               teststat = "score",
                               score.cutoff = params$score.cutoff,
                               use.SPA.score = FALSE,
                               use.marg.p = use.marg.p,
                               alpha=alpha)
    
    res_score_wSPA_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                                    Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                    sim_map=sim_map,
                                    glm_input = attrDat$nulldata,
                                    teststat = "score",
                                    score.cutoff = params$score.cutoff,
                                    use.SPA.score = TRUE,
                                    use.marg.p = use.marg.p,
                                    alpha=alpha)
    
    
    #Run algorithm with method DD and FET statistic
    res_FET =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                           Gmat_ctrl=attrDat$geno_ctrl_matrix,
                           sim_map=sim_map,
                           glm_input = attrDat$nulldata,
                           teststat = "FET",
                           use.marg.p = use.marg.p,
                           L=L,
                           alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                             Gmat_ctrl=attrDat$geno_ctrl_matrix,
                             sim_map=sim_map,
                             glm_input = attrDat$nulldata,
                             teststat = "score",
                             score.cutoff = params$score.cutoff,
                             use.SPA.score = FALSE,
                             use.marg.p = use.marg.p,
                             L=L,
                             alpha=alpha)
    
    res_score_wSPA =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                                  Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                  sim_map=sim_map,
                                  glm_input = attrDat$nulldata,
                                  teststat = "score",
                                  score.cutoff = params$score.cutoff,
                                  use.SPA.score = TRUE,
                                  use.marg.p = use.marg.p,
                                  L=L,
                                  alpha=alpha)
    
    #Performance
    
    #Ground truth (leaves)
    m = sim_map$num_snps
    nonnulls = sim_map$alt_snps
    nonnulls_omni = sim_map$alt_omni_snps
    nulls = sim_map$null_snps
    
    #Ground truth - domain collapsing
    m_dc = length(sim_map$domain_snp_list)
    nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
    nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
    
    for (l in 1:L){
      
      #Method FET + DD
      Sim.mat[1,i,1,l,1] = sum(res_FET$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,1,l,2] = sum(res_FET$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,1,l,3] = sum(res_FET$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,1,l,4] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,1,l,5] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,1,l,6] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls_omni) #FN omni
      Sim.mat[1,i,1,l,7] = res_FET$LayerSumm[[l]]$FDP.l
      
      #Method score + noSPA + DD
      Sim.mat[1,i,2,l,1] = sum(res_score$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,2,l,2] = sum(res_score$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,2,l,3] = sum(res_score$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,2,l,4] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,2,l,5] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,2,l,6] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls_omni) #FN omni
      Sim.mat[1,i,2,l,7] = res_score$LayerSumm[[l]]$FDP.l
      
      #Method scoreSPA + DD
      Sim.mat[1,i,3,l,1] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,3,l,2] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,3,l,3] = sum(res_score_wSPA$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,3,l,4] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,3,l,5] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,3,l,6] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls_omni) #FN omni
      Sim.mat[1,i,3,l,7] = res_score_wSPA$LayerSumm[[l]]$FDP.l
      
      #Method FET + DD + domain collapsing (Power comparison only)
      Sim.mat[1,i,4,l,1] = sum(res_FET$S.list.dc[[l]]%in%nonnulls_dc) #TP
      Sim.mat[1,i,4,l,2] = sum(res_FET$S.list.dc[[l]]%in%nulls_dc) #FP
      Sim.mat[1,i,4,l,3] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nulls_dc) #TN
      Sim.mat[1,i,4,l,4] = sum(setdiff(seq(m_dc),res_FET$S.list.dc[[l]])%in%nonnulls_dc) #FN
      
      #Method score + DD + domain collapsing (Power comparison only)
      Sim.mat[1,i,5,l,1] = sum(res_score$S.list.dc[[l]]%in%nonnulls_dc) #TP
      Sim.mat[1,i,5,l,2] = sum(res_score$S.list.dc[[l]]%in%nulls_dc) #FP
      Sim.mat[1,i,5,l,3] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nulls_dc) #TN
      Sim.mat[1,i,5,l,4] = sum(setdiff(seq(m_dc),res_score$S.list.dc[[l]])%in%nonnulls_dc) #FN
      
      #Method score + DD + domain collapsing (Power comparison only)
      Sim.mat[1,i,6,l,1] = sum(res_score_wSPA$S.list.dc[[l]]%in%nonnulls_dc) #TP
      Sim.mat[1,i,6,l,2] = sum(res_score_wSPA$S.list.dc[[l]]%in%nulls_dc) #FP
      Sim.mat[1,i,6,l,3] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nulls_dc) #TN
      Sim.mat[1,i,6,l,4] = sum(setdiff(seq(m_dc),res_score_wSPA$S.list.dc[[l]])%in%nonnulls_dc) #FN
    }
    
    #Method FET + domain collapsing
    Sim.mat[1,i,7,1,1] = sum(res_FET_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,7,1,2] = sum(res_FET_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,7,1,3] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,7,1,4] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nonnulls_dc) #FN
    
    #Method score + domain collapsing
    Sim.mat[1,i,8,1,1] = sum(res_score_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,8,1,2] = sum(res_score_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,8,1,3] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,8,1,4] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nonnulls_dc) #FN
    
    #Method score + domain collapsing
    Sim.mat[1,i,9,1,1] = sum(res_score_wSPA_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,9,1,2] = sum(res_score_wSPA_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,9,1,3] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,9,1,4] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nonnulls_dc) #FN
    
  }
  
  return(Sim.mat)
}


run_DGP2_sing2 <- function(N1,
                          sampsize.factors,
                          sim_map,
                          params,
                          midp=TRUE,
                          use.marg.p=FALSE,
                          covars=FALSE,
                          alpha,L,seed){
  
  set.seed(seed)
  
  Sim.mat = array(NA,c(1,length(sampsize.factors),6,L,6)) #1 iteration, prevalence, 3 teststats x 2 collapsing strategies, L layers, and 6 performance metrics
  
  
  ### Generate Sample Data ###
  p0 = 1/100
  
  for(i in seq_along(sampsize.factors)){
    
    ## Generate case-control samples
    N1 <- N1*sampsize.factors[i]
    N0 <- N1
    
    ## Create data
    attrDat <- DGP2(Npop=ceiling(3*N1/p0),
                    N0 = N0,N1 = N1,
                    sim_map = sim_map,
                    covars=covars,
                    params = params,
                    p0 = p0)
    
    res_FET_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                             Gmat_ctrl=attrDat$geno_ctrl_matrix,
                             sim_map=sim_map,
                             glm_input = attrDat$nulldata,
                             teststat = "FET",
                             use.marg.p = use.marg.p,
                             alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                               Gmat_ctrl=attrDat$geno_ctrl_matrix,
                               sim_map=sim_map,
                               glm_input = attrDat$nulldata,
                               teststat = "score",
                               score.cutoff = params$score.cutoff,
                               use.SPA.score = FALSE,
                               use.marg.p = use.marg.p,
                               alpha=alpha)
    
    res_score_wSPA_dc =SATET_sim_dc(Gmat_case=attrDat$geno_case_matrix,
                                    Gmat_ctrl=attrDat$geno_ctrl_matrix,
                                    sim_map=sim_map,
                                    glm_input = attrDat$nulldata,
                                    teststat = "score",
                                    score.cutoff = params$score.cutoff,
                                    use.SPA.score = TRUE,
                                    use.marg.p = use.marg.p,
                                    alpha=alpha)
    
    
    #Run algorithm with method DD and FET statistic
    res_FET =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                       Gmat_ctrl=attrDat$geno_ctrl_matrix,
                       agg_snps = FALSE,
                       sim_map=sim_map,
                       glm_input = attrDat$nulldata,
                       teststat = "FET",
                       use.marg.p = use.marg.p,
                       L=L,
                       alpha=alpha)
    
    #Run algorithm with method DD and score statistic (logistic regression)
    res_score =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                         Gmat_ctrl=attrDat$geno_ctrl_matrix,
                         agg_snps = FALSE,
                         sim_map=sim_map,
                         glm_input = attrDat$nulldata,
                         teststat = "score",
                         score.cutoff = params$score.cutoff,
                         use.SPA.score = FALSE,
                         use.marg.p = use.marg.p,
                         L=L,
                         alpha=alpha)
    
    res_score_wSPA =SATET_sim_snp(Gmat_case=attrDat$geno_case_matrix,
                              Gmat_ctrl=attrDat$geno_ctrl_matrix,
                              agg_snps = FALSE,
                              sim_map=sim_map,
                              glm_input = attrDat$nulldata,
                              teststat = "score",
                              score.cutoff = params$score.cutoff,
                              use.SPA.score = TRUE,
                              use.marg.p = use.marg.p,
                              L=L,
                              alpha=alpha)
    
    #Performance
    
    #Ground truth (leaves)
    m = sim_map$num_snps
    nonnulls = sim_map$alt_snps
    nonnulls_omni = sim_map$alt_omni_snps
    nulls = sim_map$null_snps
    
    #Ground truth - domain collapsing
    m_dc = length(sim_map$domain_snp_list)
    nonnulls_dc = unlist(sim_map$domain_w_sa_alt_snps)
    nulls_dc = setdiff(seq(m_dc),nonnulls_dc)
    
    for (l in 1:L){
      
      #Method FET + DD
      Sim.mat[1,i,1,l,1] = sum(res_FET$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,1,l,2] = sum(res_FET$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,1,l,3] = sum(res_FET$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,1,l,4] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,1,l,5] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,1,l,6] = sum(setdiff(seq(m),res_FET$S.list[[l]])%in%nonnulls_omni) #FN omni
      
      #Method score + noSPA + DD
      Sim.mat[1,i,2,l,1] = sum(res_score$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,2,l,2] = sum(res_score$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,2,l,3] = sum(res_score$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,2,l,4] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,2,l,5] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,2,l,6] = sum(setdiff(seq(m),res_score$S.list[[l]])%in%nonnulls_omni) #FN omni
      
      #Method scoreSPA + DD
      Sim.mat[1,i,3,l,1] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls) #TP
      Sim.mat[1,i,3,l,2] = sum(res_score_wSPA$S.list[[l]]%in%nonnulls_omni) #TP omni
      Sim.mat[1,i,3,l,3] = sum(res_score_wSPA$S.list[[l]]%in%nulls) #FP
      Sim.mat[1,i,3,l,4] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nulls) #TN
      Sim.mat[1,i,3,l,5] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls) #FN
      Sim.mat[1,i,3,l,6] = sum(setdiff(seq(m),res_score_wSPA$S.list[[l]])%in%nonnulls_omni) #FN omni
      
    }
    
    #Method FET + DD + domain collapsing
    Sim.mat[1,i,4,1,1] = sum(res_FET_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,4,1,2] = sum(res_FET_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,4,1,3] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,4,1,4] = sum(setdiff(seq(m_dc),res_FET_dc$S.list.dc)%in%nonnulls_dc) #FN
    
    #Method score + DD + domain collapsing
    Sim.mat[1,i,5,1,1] = sum(res_score_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,5,1,2] = sum(res_score_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,5,1,3] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,5,1,4] = sum(setdiff(seq(m_dc),res_score_dc$S.list.dc)%in%nonnulls_dc) #FN
    
    #Method score + DD + domain collapsing
    Sim.mat[1,i,6,1,1] = sum(res_score_wSPA_dc$S.list.dc%in%nonnulls_dc) #TP
    Sim.mat[1,i,6,1,2] = sum(res_score_wSPA_dc$S.list.dc%in%nulls_dc) #FP
    Sim.mat[1,i,6,1,3] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nulls_dc) #TN
    Sim.mat[1,i,6,1,4] = sum(setdiff(seq(m_dc),res_score_wSPA_dc$S.list.dc)%in%nonnulls_dc) #FN
    
    
  }
  
  return(Sim.mat)
}


simGenomicStructure_global_omnigenic_random <- function(snp_dat,
                                                        snp_clust_len,
                                                        num_rep_clust){
  require(data.table)
  
  #Get random SNPset of contiguous SNPs
  getrandom_contig_snpset <- function(v,p){
    possibleIndex = seq(length(v) - p + 1)
    firstIndex = sample(possibleIndex, 1)
    out = v[firstIndex:(firstIndex + p -1)]
    #Check for contiguous SNPs
    while(any(rle(diff(out))$values>1)){
      firstIndex = sample(possibleIndex, 1)
      out = v[firstIndex:(firstIndex + p -1)]
    }
    return(out)
  }
  
  #Renumber variants to use less memory
  domain_snp_list = split(snp_dat$loc_adj,#rank_loc_adj,
                          as.character(snp_dat$subRVIS.Domain.Name))
  
  domain_unique_snp_list = unname(lapply(domain_snp_list,unique))  
  
  #Get snp IDs for each cluster size
  domain_w_sa_alt_snps <-vector("list",length(snp_clust_len))
  snps_in_cluster <- vector("list",length(snp_clust_len))
  
  omni_snps_in_cluster <- vector("list",length(snp_clust_len))
  
  #For each domain keep track of the snp indices - immutable
  domain_snp_list <- relist(seq_along(unlist(domain_unique_snp_list)),
                            skeleton = domain_unique_snp_list)
  
  domain_snp_list_ <- domain_snp_list
  
  for(k in rev(seq_along(snp_clust_len))){
    
    #First randomly select domains for strong association SNPs
    domain_sa_snp_len_k <- sample(which(unname(lapply(domain_snp_list_,
                                                      function(x) length(x))>=5*snp_clust_len[k])),
                                  num_rep_clust[k])
    
    
    snps_in_cluster[[k]] <- unlist(lapply(domain_snp_list_[domain_sa_snp_len_k],
                                          function(x) getrandom_contig_snpset(x,snp_clust_len[k])))
    
    #Update domain/snp info
    domain_snp_list_ <- lapply(domain_snp_list,function(x) setdiff(x,unlist(snps_in_cluster)))
    
    domain_wa_snp_len_k <- sample(which(unname(lapply(domain_snp_list_,
                                                      function(x) length(x))>=2*snp_clust_len[k])),
                                  num_rep_clust[k])
    
    omni_snps_in_cluster[[k]] <- unlist(lapply(domain_snp_list_[domain_wa_snp_len_k],
                                               function(x) getrandom_contig_snpset(x,snp_clust_len[k])))
    
    domain_snp_list_ <- lapply(domain_snp_list,function(x) setdiff(x, unlist(omni_snps_in_cluster)))
    
    domain_w_sa_alt_snps[[k]] <- domain_sa_snp_len_k
    
  }
  
  alt_snps = unlist(snps_in_cluster)
  alt_omni_snps = unlist(omni_snps_in_cluster)
  null_snps = setdiff(unlist(domain_snp_list),
                      c(alt_snps,alt_omni_snps))
  
  return(list("type"="global",
              "num_snps"=length(unique(snp_dat$loc_adj)),
              "domain_snp_list"=domain_snp_list,
              "domain_w_sa_alt_snps"=domain_w_sa_alt_snps,
              "alt_snps"=alt_snps,
              "alt_omni_snps"=alt_omni_snps,
              "null_snps"=null_snps,
              "snp_clust_len"=snp_clust_len,
              "num_rep_clust"=num_rep_clust,
              "snps_in_cluster"=snps_in_cluster,
              "omni_snps_in_cluster"=omni_snps_in_cluster))
}

simGenomicStructure_global_omnigenic_random2 <- function(snp_dat,
                                                        snp_clust_len,
                                                        num_rep_clust,
                                                        omni_prop){
  require(data.table)
  
  #Get random SNPset of contiguous SNPs
  getrandom_contig_snpset <- function(v,p){
    possibleIndex = seq(length(v) - p + 1)
    firstIndex = sample(possibleIndex, 1)
    out = v[firstIndex:(firstIndex + p -1)]
    #Check for contiguous SNPs
    while(any(rle(diff(out))$values>1)){
      firstIndex = sample(possibleIndex, 1)
      out = v[firstIndex:(firstIndex + p -1)]
    }
    return(out)
  }
  
  #Renumber variants to use less memory
  domain_snp_list = split(snp_dat$loc_adj,#rank_loc_adj,
                          as.character(snp_dat$subRVIS.Domain.Name))
  
  domain_unique_snp_list = unname(lapply(domain_snp_list,unique))  
  
  #Get snp IDs for each cluster size
  domain_w_sa_alt_snps <-vector("list",length(snp_clust_len))
  snps_in_cluster <- vector("list",length(snp_clust_len))
  
  omni_snps_in_cluster <- vector("list",length(snp_clust_len))
  
  #For each domain keep track of the snp indices - immutable
  domain_snp_list <- relist(seq_along(unlist(domain_unique_snp_list)),
                            skeleton = domain_unique_snp_list)
  
  domain_snp_list_ <- domain_snp_list
  
  for(k in rev(seq_along(snp_clust_len))){
    
    #First randomly select domains for strong association SNPs
    domain_sa_snp_len_k <- sample(which(unname(lapply(domain_snp_list_,
                                                      function(x) length(x))>=5*snp_clust_len[k])),
                                  num_rep_clust[k])
    
    
    snps_in_cluster[[k]] <- unlist(lapply(domain_snp_list_[domain_sa_snp_len_k],
                                          function(x) getrandom_contig_snpset(x,snp_clust_len[k])))
    
    #Update domain/snp info
    domain_w_sa_alt_snps[[k]] <- domain_sa_snp_len_k
    
  }
  
  # omni_snps_in_cluster[[k]] <- unlist(lapply(setdiff(domain_snp_list[domain_snp_len_k],
  #                                      snps_in_cluster[[k]]),
  #                                     function(x) sample(x,ceiling(length(x)*omni_prop),replace=FALSE)))
  
  
  
  alt_snps = unlist(snps_in_cluster)
  remaining_snps = setdiff(unlist(domain_snp_list),alt_snps)
  alt_omni_snps = sample(remaining_snps, ceiling(length(remaining_snps)*omni_prop),replace = FALSE)
  null_snps = setdiff(unlist(domain_snp_list),
                      c(alt_snps,alt_omni_snps))
  
  return(list("type"="global",
              "num_snps"=length(unique(snp_dat$loc_adj)),
              "domain_snp_list"=domain_snp_list,
              "domain_w_sa_alt_snps"=domain_w_sa_alt_snps,
              "alt_snps"=alt_snps,
              "alt_omni_snps"=alt_omni_snps,
              "null_snps"=null_snps,
              "snp_clust_len"=snp_clust_len,
              "num_rep_clust"=num_rep_clust,
              "snps_in_cluster"=snps_in_cluster))
}



simGenomicStructure_local_omnigenic_random <- function(snp_dat,
                                                       snp_clust_len,
                                                       num_rep_clust){
  require(data.table)
  
  #set.seed(seed) #only used for sampling from domain and variants - can be changed
  
  #From snp_dat, create list of UNIQUE variant vectors for each domain
  #Each element of the list is a domain with the corresponding genomic positions found in the data
  
  #Renumber variants to use less memory
  domain_snp_list = split(snp_dat$loc_adj,#rank_loc_adj,
                          as.character(snp_dat$subRVIS.Domain.Name))
  
  domain_unique_snp_list = unname(lapply(domain_snp_list,unique))  
  
  #Get snp IDs for each cluster size
  domain_w_sa_alt_snps <-vector("list",length(snp_clust_len))
  snps_in_cluster <- vector("list",length(snp_clust_len))
  
  omni_snps_in_cluster <- vector("list",length(snp_clust_len))
  
  #For each domain keep track of the snp indices
  domain_snp_list <- relist(seq_along(unlist(domain_unique_snp_list)),
                            skeleton = domain_unique_snp_list)
  
  for(k in rev(seq_along(snp_clust_len))){
    #Ensure that unique domains are selected
    domain_snp_len_k <- sample(which(unname(lapply(domain_unique_snp_list,
                                                   function(x) length(x))>=5*snp_clust_len[k])),
                               num_rep_clust[k])
    
    while(k<length(snp_clust_len) & 
          any(domain_snp_len_k%in%unlist(domain_w_sa_alt_snps[(k+1):length(snp_clust_len)]))){
      #Resample
      domain_snp_len_k <- sample(which(unname(lapply(domain_unique_snp_list,
                                                     function(x) length(x))>=5*snp_clust_len[k])),
                                 num_rep_clust[k])
    }
    
    domain_w_sa_alt_snps[[k]] <- domain_snp_len_k
    
    #Identify all snps within each selected domain in domain_snp_len_k
    getrandom_snpset <- function(v,p){
      possibleIndex = seq(length(v) - p + 1)
      firstIndex = sample(possibleIndex, 1)
      v[firstIndex:(firstIndex + p -1)]
    }
    
    snps_in_cluster[[k]] <- unlist(lapply(domain_snp_list[domain_w_sa_alt_snps[[k]]],
                                          function(x) getrandom_snpset(x,snp_clust_len[k])))
    
    #For omnigenic variants - within selected domains with non-omnigenic snpsets,
    #set all other snps to alternatives
    #store variants outside the non-omnigenic snp sets
    
    omni_snps_in_cluster[[k]] <- setdiff(unlist(domain_snp_list[domain_snp_len_k]),
                                         snps_in_cluster[[k]])
    
  }
  
  alt_snps = unlist(snps_in_cluster)
  alt_omni_snps = unlist(omni_snps_in_cluster)
  null_snps = setdiff(unlist(domain_snp_list),
                      c(alt_snps,alt_omni_snps))
  
  return(list("type"="local",
              "num_snps"=length(unique(snp_dat$loc_adj)),
              "domain_snp_list"=domain_snp_list,
              "domain_w_sa_alt_snps"=domain_w_sa_alt_snps,
              "alt_snps"=alt_snps,
              "alt_omni_snps"=alt_omni_snps,
              "null_snps"=null_snps,
              "snp_clust_len"=snp_clust_len,
              "num_rep_clust"=num_rep_clust,
              "snps_in_cluster"=snps_in_cluster,
              "omni_snps_in_cluster"=omni_snps_in_cluster))
}

#Generate positions with varying omnigenic proportions
simGenomicStructure_local_omnigenic_random2 <- function(snp_dat,
                                                       snp_clust_len,
                                                       omni_prop,
                                                       num_rep_clust){
  require(data.table)

  #set.seed(seed) #only used for sampling from domain and variants - can be changed

  #From snp_dat, create list of UNIQUE variant vectors for each domain
  #Each element of the list is a domain with the corresponding genomic positions found in the data

  #Renumber variants to use less memory
  domain_snp_list = split(snp_dat$loc_adj,#rank_loc_adj,
                          as.character(snp_dat$subRVIS.Domain.Name))

  domain_unique_snp_list = unname(lapply(domain_snp_list,unique))

  #Get snp IDs for each cluster size
  domain_w_sa_alt_snps <-vector("list",length(snp_clust_len))
  snps_in_cluster <- vector("list",length(snp_clust_len))

  omni_snps_in_cluster <- vector("list",length(snp_clust_len))

  #For each domain keep track of the snp indices
  domain_snp_list <- relist(seq_along(unlist(domain_unique_snp_list)),
                            skeleton = domain_unique_snp_list)

  for(k in rev(seq_along(snp_clust_len))){
    #Ensure that unique domains are selected
    domain_snp_len_k <- sample(which(unname(lapply(domain_unique_snp_list,
                                                   function(x) length(x))>=5*snp_clust_len[k])),
                               num_rep_clust[k])

    while(k<length(snp_clust_len) &
          any(domain_snp_len_k%in%unlist(domain_w_sa_alt_snps[(k+1):length(snp_clust_len)]))){
      #Resample
      domain_snp_len_k <- sample(which(unname(lapply(domain_unique_snp_list,
                                                     function(x) length(x))>=5*snp_clust_len[k])),
                                 num_rep_clust[k])
    }

    domain_w_sa_alt_snps[[k]] <- domain_snp_len_k

    #Identify all snps within each selected domain in domain_snp_len_k
    getrandom_snpset <- function(v,p){
      possibleIndex = seq(length(v) - p + 1)
      firstIndex = sample(possibleIndex, 1)
      v[firstIndex:(firstIndex + p -1)]
    }

    snps_in_cluster[[k]] <- unlist(lapply(domain_snp_list[domain_w_sa_alt_snps[[k]]],
                                          function(x) getrandom_snpset(x,snp_clust_len[k])))

    #For omnigenic variants - within selected domains with non-omnigenic snpsets,
    #set all other snps to alternatives
    #store variants outside the non-omnigenic snp sets

    remaining_snps_in_cluster_k <- setdiff(unlist(domain_snp_list[domain_snp_len_k]),
                                           snps_in_cluster[[k]])
    
    omni_snps_in_cluster[[k]] <- sample(remaining_snps_in_cluster_k,
                                        ceiling(length(remaining_snps_in_cluster_k)*omni_prop),
                                        replace = FALSE)

  }

  
  alt_snps = unlist(snps_in_cluster)
  alt_omni_snps = unlist(omni_snps_in_cluster)
  null_snps = setdiff(unlist(domain_snp_list),
                      c(alt_snps,alt_omni_snps))
  
  return(list("type"="local",
              "num_snps"=length(unique(snp_dat$loc_adj)),
              "domain_snp_list"=domain_snp_list,
              "domain_w_sa_alt_snps"=domain_w_sa_alt_snps,
              "alt_snps"=alt_snps,
              "alt_omni_snps"=alt_omni_snps,
              "null_snps"=null_snps,
              "snp_clust_len"=snp_clust_len,
              "num_rep_clust"=num_rep_clust,
              "snps_in_cluster"=snps_in_cluster,
              "omni_snps_in_cluster"=omni_snps_in_cluster))
  
}

genSpMat <- function(nrows, ncols, col_probs) {
  require(Matrix)
  r <- lapply(1:ncols, function(x) {
    p <- col_probs[x]
    i <- sample.int(2L, size = nrows, replace = T, prob = c(1 - p, p))
    which(i == 2L)
  })
  rl <- lengths(r)
  nc <- rep(1:ncols, times = rl) # col indexes
  nr <- unlist(r) # row index
  ddims <- c(nrows, ncols)
  #ngTMatrix output
  sparseMatrix(i = nr, j = nc, dims = ddims,giveCsparse = FALSE)
}

DGP1 <- function(Npop,
                 N0=500,
                 N1=500,
                 sim_map,
                 covars=FALSE,
                 params,
                 p0=0.01){

  require(Matrix)

  n_vars <- length(unlist(sim_map$domain_snp_list))

  #Generate covariates
  X1 <- rnorm(Npop,0,1)
  X2 <- rbinom(Npop,1,0.5)

  #Generate sparse genotype matrix
  #Gmat <- genSpMat2(Npop,n_vars,params$maf)
  Gmat <- genSpMat(Npop,n_vars,params$maf)

  #Alternative SNP locations
  altSNPloc <- c(sim_map$alt_snps,sim_map$alt_omni_snps)

  snp_clust_len <- sim_map$snp_clust_len
  num_rep_clust <- sim_map$num_rep_clust
  num_rep_omni_clust <- sim_map$num_rep_omni_clust

  betas = rep(0,n_vars)

  for(i in seq_along(snp_clust_len)){
      betas[sim_map$snps_in_cluster[[i]]] <- params$betas_alt[[i]]/params$factor.alt[i]
      betas[sim_map$omni_snps_in_cluster[[i]]] <- params$betas_alt_omni[[i]]/params$factor.alt.omni[i]
  }

  #Starting value for alpha0 
  alpha0 <- log(p0/(1-p0)) 
  
  pr <- rep(1,Npop)
  #Covariates and G'beta
  Gb <- as.vector(Gmat[,altSNPloc] %*% betas[altSNPloc])
  Xg_Gb <- params$gamma[1]*X1 + params$gamma[2]*X2 + Gb
  
  #Starting value for pr
  pr <- plogis(alpha0 + Xg_Gb)
  
  #Adjust intercept until population prevalence is reached (to a certain tolerance)
  if(covars){
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  } else{
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  }
  
  Y <- rbinom(Npop,1,pr)

  #Randomly select N1 cases and N0 controls
  #ngTMatrix objects
  sampcase <- sample(which(Y==1),N1)
  sampctrl <- sample(which(Y==0),N0)

  #Get attribute matrices
  geno_case_matrix <- Gmat[sampcase,]
  geno_ctrl_matrix <- Gmat[sampctrl,]

  X1ctrl <- X1[sampctrl]
  X1case <- X1[sampcase]
  X2ctrl <- X2[sampctrl]
  X2case <- X2[sampcase]

  rm(Gmat)

  #Output
  if(covars){
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0)),
                             X1=c(X1case,X1ctrl),
                             X2=c(X2case,X2ctrl)),#Note order,
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  } else{
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0))),
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  }
}


DGP1_conf <- function(Npop,
                 N0=500,
                 N1=500,
                 sim_map,
                 covars=FALSE,
                 params,
                 p0=0.01){
  
  require(Matrix)
  
  n_vars <- length(unlist(sim_map$domain_snp_list))
  
  #Generate covariates
  X1 <- rnorm(Npop,0,1)
  #X2 <- rbinom(Npop,1,0.5)
  
  #Generate sparse genotype matrix
  #Gmat <- genSpMat2(Npop,n_vars,params$maf)
  Gmat <- genSpMat(Npop,n_vars,params$maf)

  #Alternative SNP locations
  altSNPloc <- c(sim_map$alt_snps,sim_map$alt_omni_snps)
  
  ##Generate Binary Confounder
  
  n_alt_samp <- 80
  n_nul_samp <- 320
  #1) sample n_alt_samp from altSNPloc and n_nul_samp from remainder  
  alt_samp <- sample(altSNPloc,n_alt_samp)
  nul_samp <- sample(setdiff(seq(n_vars),altSNPloc),n_nul_samp)
  
  #2) For each individual, count number of locations linked to at least one rare allele
  A_samp_i <- rowSums(Gmat[,c(alt_samp,nul_samp)])
  quant_A_samp <- quantile(A_samp_i,0.8)
  
  #3) For each individual i, if A_samp_i > quant_A_samp, set X2 = 1; else X2=0
  X2 <- ifelse(A_samp_i>quant_A_samp,1,0)
  
  snp_clust_len <- sim_map$snp_clust_len
  num_rep_clust <- sim_map$num_rep_clust
  num_rep_omni_clust <- sim_map$num_rep_omni_clust
  
  betas = rep(0,n_vars)
  
  for(i in seq_along(snp_clust_len)){
    betas[sim_map$snps_in_cluster[[i]]] <- params$betas_alt[[i]]/params$factor.alt[i]
    betas[sim_map$omni_snps_in_cluster[[i]]] <- params$betas_alt_omni[[i]]/params$factor.alt.omni[i]
  }
  
  #Starting value for alpha0 
  alpha0 <- log(p0/(1-p0)) 
  
  pr <- rep(1,Npop)
  #Covariates and G'beta
  Gb <- as.vector(Gmat[,altSNPloc] %*% betas[altSNPloc])
  Xg_Gb <- params$gamma[1]*X1 + params$gamma[2]*X2 + Gb
  
  #Starting value for pr
  pr <- plogis(alpha0 + Xg_Gb)
  
  #Adjust intercept until population prevalence is reached (to a certain tolerance)
  if(covars){
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  } else{
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  }
  
  Y <- rbinom(Npop,1,pr)
  
  #Randomly select N1 cases and N0 controls
  #ngTMatrix objects
  sampcase <- sample(which(Y==1),N1)
  sampctrl <- sample(which(Y==0),N0)
  
  #Get attribute matrices
  geno_case_matrix <- Gmat[sampcase,]
  geno_ctrl_matrix <- Gmat[sampctrl,]
  
  X1ctrl <- X1[sampctrl]
  X1case <- X1[sampcase]
  X2ctrl <- X2[sampctrl]
  X2case <- X2[sampcase]
  
  rm(Gmat)
  
  #Output
  if(covars){
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0)),
                             X1=c(X1case,X1ctrl),
                             X2=c(X2case,X2ctrl)),#Note order,
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  } else{
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0))),
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  }
}


DGP2 <- function(Npop,
                 N0=500, 
                 N1=500, 
                 sim_map,
                 covars=FALSE,
                 params, 
                 p0=0.01){
  
  require(Matrix)
  
  
  n_vars <- length(unlist(sim_map$domain_snp_list))
  
  #Generate covariates
  X1 <- rnorm(Npop,0,1)
  X2 <- rbinom(Npop,1,0.5)
  
  #Generate sparse genotype matrix
  Gmat <- genSpMat(Npop,n_vars,params$maf)
  
  #Alternative SNP locations
  altSNPloc <- sim_map$alt_snps
  alt_omniSNPloc <- sim_map$alt_omni_snps
  
  #Generate betas's for causal, non-omnigenic variants
  betas <- rep(0,n_vars)
  
  betas[altSNPloc] <- log(-params$beta_rho_factor[1]*log10(params$maf[altSNPloc]))
  betas[alt_omniSNPloc] <- log(-params$beta_rho_factor[2]*log10(params$maf[alt_omniSNPloc]))
  
  #Create logistic model of phenotype
  # mean_Xg_Gb = median(params$gamma[1]*X1 + params$gamma[2]*X2 + 
  #                  as.vector(Gmat[,c(altSNPloc,alt_omniSNPloc)] %*% betas[c(altSNPloc,alt_omniSNPloc)]))
  
  #Starting value for alpha0 
  alpha0 <- log(p0/(1-p0)) 
  
  #Covariates and G'beta
  Gb <- as.vector(Gmat[,c(altSNPloc,alt_omniSNPloc)] %*% betas[c(altSNPloc,alt_omniSNPloc)])
  Xg_Gb <- params$gamma[1]*X1 + params$gamma[2]*X2 + Gb
  
  #Starting value for pr
  pr <- plogis(alpha0 + Xg_Gb)
  
  #Adjust intercept until population prevalence is reached (to a certain tolerance)
  if(covars){
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  } else{
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  }

  Y <- rbinom(Npop,1,pr)
  
  #Randomly select N1 cases and N0 controls
  #ngTMatrix objects
  sampcase <- sample(which(Y==1),N1)
  sampctrl <- sample(which(Y==0),N0)
  
  #Get attribute matrices
  geno_case_matrix <- Gmat[sampcase,]
  geno_ctrl_matrix <- Gmat[sampctrl,]
  
  X1ctrl <- X1[sampctrl]
  X1case <- X1[sampcase]
  X2ctrl <- X2[sampctrl]
  X2case <- X2[sampcase]
  
  rm(Gmat)
  
  #Output
  if(covars){
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0)),
                             X1=c(X1case,X1ctrl),
                             X2=c(X2case,X2ctrl)),#Note order,
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  } else{
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0))),
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  }
  
}

DGP2_conf <- function(Npop,
                 N0=500, 
                 N1=500, 
                 sim_map,
                 covars=FALSE,
                 params, 
                 p0=0.01){
  
  require(Matrix)
  
  
  n_vars <- length(unlist(sim_map$domain_snp_list))
  
  #Generate covariates
  X1 <- rnorm(Npop,0,1)
  #X2 <- rbinom(Npop,1,0.5)
  
  #Generate sparse genotype matrix
  Gmat <- genSpMat(Npop,n_vars,params$maf)
  
  #Alternative SNP locations
  altSNPloc <- sim_map$alt_snps
  alt_omniSNPloc <- sim_map$alt_omni_snps
  
  ##Generate Binary Confounder
  
  n_alt_samp <- 80
  n_nul_samp <- 320
  #1) sample n_alt_samp from altSNPloc and n_nul_samp from remainder  
  alt_samp <- sample(altSNPloc,n_alt_samp)
  nul_samp <- sample(setdiff(seq(n_vars),altSNPloc),n_nul_samp)
  
  #2) For each individual, count number of locations linked to at least one rare allele
  A_samp_i <- rowSums(Gmat[,c(alt_samp,nul_samp)])
  quant_A_samp <- quantile(A_samp_i,0.8)
  
  #3) For each individual i, if A_samp_i > quant_A_samp, set X2 = 1; else X2=0
  X2 <- ifelse(A_samp_i>quant_A_samp,1,0)
  
  #Generate betas's for causal, non-omnigenic variants
  betas <- rep(0,n_vars)
  
  betas[altSNPloc] <- log(-params$beta_rho_factor[1]*log10(params$maf[altSNPloc]))
  betas[alt_omniSNPloc] <- log(-params$beta_rho_factor[2]*log10(params$maf[alt_omniSNPloc]))
  
  #Create logistic model of phenotype
  # mean_Xg_Gb = median(params$gamma[1]*X1 + params$gamma[2]*X2 + 
  #                  as.vector(Gmat[,c(altSNPloc,alt_omniSNPloc)] %*% betas[c(altSNPloc,alt_omniSNPloc)]))
  
  #Starting value for alpha0 
  alpha0 <- log(p0/(1-p0)) 
  
  #Covariates and G'beta
  Gb <- as.vector(Gmat[,c(altSNPloc,alt_omniSNPloc)] %*% betas[c(altSNPloc,alt_omniSNPloc)])
  Xg_Gb <- params$gamma[1]*X1 + params$gamma[2]*X2 + Gb
  
  #Starting value for pr
  pr <- plogis(alpha0 + Xg_Gb)
  
  #Adjust intercept until population prevalence is reached (to a certain tolerance)
  if(covars){
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Xg_Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  } else{
    if(mean(pr)>p0){
      #Decrement alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 - 0.05
      }
    } else{
      #Increment alpha0
      while(abs(mean(pr)-p0) > 10^(round(log10(p0))-1)){
        pr <- plogis(alpha0 + Gb)
        alpha0 <- alpha0 + 0.05
      }
    }
  }
  
  Y <- rbinom(Npop,1,pr)
  
  #Randomly select N1 cases and N0 controls
  #ngTMatrix objects
  sampcase <- sample(which(Y==1),N1)
  sampctrl <- sample(which(Y==0),N0)
  
  #Get attribute matrices
  geno_case_matrix <- Gmat[sampcase,]
  geno_ctrl_matrix <- Gmat[sampctrl,]
  
  X1ctrl <- X1[sampctrl]
  X1case <- X1[sampcase]
  X2ctrl <- X2[sampctrl]
  X2case <- X2[sampcase]
  
  rm(Gmat)
  
  #Output
  if(covars){
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0)),
                             X1=c(X1case,X1ctrl),
                             X2=c(X2case,X2ctrl)),#Note order,
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  } else{
    list(nulldata=data.frame(phenotype=rep(c(1,0),times=c(N1,N0))),
         geno_case_matrix=geno_case_matrix,
         geno_ctrl_matrix=geno_ctrl_matrix)
  }
  
}



FDR_summ_dc <- function(res,null_indx,alt_indx){
  ## # #print out empirical FDR for each layer
    print(paste("FDR:",
                sum(res$S.list.dc[[1]]%in%null_indx)/max(1,length(res$S.list.dc[[1]])),
                sep=" "))
    print(paste("Sensitivity:",
                sum(res$S.list.dc[[1]]%in%alt_indx)/length(alt_indx),
                sep=" "))
}


FDR_summ_omni <- function(res,null_indx,alt_indx,alt_omni_indx){
  ## # #print out empirical FDR for each layer
  for(l in seq_along(res$S.list)){
    print(paste("Layer",l,"FDR:",
                sum(res$S.list[[l]]%in%null_indx)/max(1,length(res$S.list[[l]])),
                sep=" "))
    print(paste("Layer",l,"Sensitivity (non-omni):",
                sum(res$S.list[[l]]%in%alt_indx)/length(alt_indx),
                sep=" "))
    print(paste("Layer",l,"Sensitivity (omni):",
                sum(res$S.list[[l]]%in%alt_omni_indx)/length(alt_omni_indx),
                sep=" "))
  }
}


FDR_summ_agg <- function(res){
  ## # #print out empirical FDR for each layer
  for(l in seq_along(res$S.list)){
    print(paste("Layer",l,"FDR:",
                sum(res$S.list[[l]]%in%res$LayerSumm[[1]]$null_leaves)/max(1,length(res$S.list[[l]])),
                sep=" "))
    print(paste("Layer",l,"Sensitivity:",
                sum(res$S.list[[l]]%in%res$LayerSumm[[1]]$alt_leaves)/length(res$LayerSumm[[1]]$alt_leaves),
                sep=" "))
  }
}

FDR_summ_sing <- function(res){
  ## # #print out empirical FDR for each layer
  for(l in seq_along(res$S.list)){
    print(paste("Layer",l,"FDR:",
                sum(res$S.list[[l]]%in%res$LayerSumm[[1]]$null_snps)/max(1,length(res$S.list[[l]])),
                sep=" "))
    print(paste("Layer",l,"Sensitivity:",
                sum(res$S.list[[l]]%in%res$LayerSumm[[1]]$alt_snps)/length(res$LayerSumm[[1]]$alt_snps),
                sep=" "))
  }
}



evaluate.sim <- function(sim.mat){
  RejAlt = apply(sim.mat[,,1],2,mean,na.rm=TRUE) #TP
  RejNul = apply(sim.mat[,,2],2,mean,na.rm=TRUE) #FP
  AccAlt = apply(sim.mat[,,4],2,mean,na.rm=TRUE) #FN
  AccNul = apply(sim.mat[,,3],2,mean,na.rm=TRUE) #TN
  Alt = apply(sim.mat[,,1]+sim.mat[,,4],2,mean,na.rm=TRUE)
  FDP = apply(sim.mat[,,2]/pmax(sim.mat[,,1]+sim.mat[,,2],1),2,mean,na.rm=TRUE)
  Sens = apply(sim.mat[,,1]/pmax(sim.mat[,,1]+sim.mat[,,4],1),2,mean,na.rm=TRUE)
  return(list("RejAlt"=RejAlt,
              "RejNul"=RejNul, "AccAlt"=AccAlt,
              "AccNul"=AccNul,"AvgSens"=Sens,
              "AvgFDP"=FDP,"NumAlt"=Alt))
}

evaluate.sim.dc <- function(sim.mat){
  RejStr = apply(sim.mat[,,1],2,mean,na.rm=TRUE) #TP
  RejNul = apply(sim.mat[,,2],2,mean,na.rm=TRUE) #FP
  AccStr = apply(sim.mat[,,4],2,mean,na.rm=TRUE) #FN
  AccNul = apply(sim.mat[,,3],2,mean,na.rm=TRUE) #TN
  AltStr = apply(sim.mat[,,1]+sim.mat[,,4],2,mean,na.rm=TRUE)
  FDP = apply(sim.mat[,,2]/pmax(sim.mat[,,1]+sim.mat[,,2],1),2,mean,na.rm=TRUE)
  SensStr = apply(sim.mat[,,1]/pmax(sim.mat[,,1]+sim.mat[,,4],1),2,mean,na.rm=TRUE)
  AltWk = NA
  RejWk = NA
  AccWk = NA
  SensWk = NA
  return(list("RejStr"=RejStr,"RejWk"=RejWk,"RejNul"=RejNul, "AccStr"=AccStr,
              "AccWk"=AccWk,"AccNul"=AccNul,
              "AvgSensStr"=SensStr,"AvgSensWk"=SensWk,
              "AvgFDP"=FDP,"NumAltStr"=AltStr,"NumAltWk"=AltWk))
}



qqPlot <- function(pval) {
  require(ggplot2)
  pval <- pval[!is.na(pval)]
  n <- length(pval)
  x <- 1:n
  dat <- data.frame(obs=sort(pval),
                    exp=x/n,
                    upper=qbeta(0.025, x, rev(x)),
                    lower=qbeta(0.975, x, rev(x)))
  
  ggplot(dat, aes(-log10(exp), -log10(obs))) +
    geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
    geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    theme_bw()
}    

