library(R.matlab)
source("./reacmat.R")

#lvm <- readMat("./Eco/Lvm_N(50)an(1.300000e-01)c(0)cp(0)mu(0)cm(0)an(0)ind(1).mat")


lvm2rn <- function(lvm,tot_org=F){
  rn <- list()
  basal <- which(lvm$tl==0)
  rn$sp.idn <- as.character(1:length(lvm$r))
  rn$sp.name <- character(length(lvm$r))
  b=1
  sp=1
  for (i in 1:length(lvm$r)){ # creation of names for the species, if the are basal o not.
    #lvm$tl represents the normalized trophic level
    if(lvm$tl[i]==0){
      #basal species
      rn$sp.name[i] <- paste("basal",b)
      b <- b+1
    }
    else{
      #non-basal species
      rn$sp.name[i] <- paste("specie",sp)
      sp <- sp+1
    }
  } 
  #rnames on variable id
  rn$sp.id <- rn$sp.name
  scod <- 1:length(rn$sp.id); names(scod) <- rn$sp.id
  
  gr2rn <- function(j){ #grow rate to reaction network
    out <- list()
    if(lvm$r[j,1]>0){ #positive grow rate
      out$reactants <- character(0)
      out["reactants"] <- list(NULL) 
      out$products <- rn$sp.name[j]
      out$r.stoichiometry <- numeric(0)
      out$p.stoichiometry <- 1
      return(out)
      } 
    else if(lvm$r[j,1]<0){ #negative grow rate
      out$reactants <- rn$sp.name[j]
      out$products <- character(0)
      out["products"] <- list(NULL)
      out$r.stoichiometry <- 1
      out$p.stoichiometry <- numeric(0)
      return(out)
    }
  }
  
  rn$reac<-lapply(which(lvm$r[,1]!=0),function(x) gr2rn(x))
  
  com2rn <- function(i,j){#community matrix interactions to reaction networks for each par (i,j)
    out <- list()
    
    if(i==j & lvm$A[i,j] != 0){#logistic par
      out$reactants <- rn$sp.name[i]
      out$products <- out$reactants
      out$r.stoichiometry <- 2
      out$p.stoichiometry <- 1
    }else if(lvm$A[i,j]>0){# comensal j over i
      out$reactants <- c(rn$sp.name[i],rn$sp.name[j])
      out$products <- out$reactants
      out$r.stoichiometry <- c(1,1)
      out$p.stoichiometry <- c(2,1)
    }else if(lvm$A[i,j]<0){#amensal j over i
      out$reactants <- c(rn$sp.name[i],rn$sp.name[j])
      out$products <- rn$sp.name[i]
      out$r.stoichiometry <- c(1,1)
      out$p.stoichiometry <- 1
    }
    return(out)
  }
  
  for(i in 1:length(lvm$r)){#geneartion of 
    print(i)
    rn$reac <- c(rn$reac,lapply(which(lvm$A[i,]!=0),function(x) com2rn(i,x)))
  }
  
  #stoiquimetric matrix generations
  mr <- mp <- matrix(0,length(lvm$r),length(rn$reac))
  rownames(mr) <- rn$sp.id; rownames(mp) <- rn$sp.id
  
  #stoiquimetric matrix assingment
  i <- 1
  for (r in rn$reac) {
    l <- length(r$reactants); lp <- length(r$products)
    j <- 1; for (s in r$reactants) { mr[s,i] <- r$r.stoichiometry[j]; j <- j + 1 } 
    j <- 1; for (s in r$products) { mp[s,i] <- r$p.stoichiometry[j]; j <- j + 1 } 
    i <- i + 1
  }
  cat(nrow(mr),"especies,",ncol(mr), "reacciones\n")
  
  if(!tot_org){
    rn <- c(rn,list(spc=scod,mr=mr,mp=mp))
    nsp <- rn.linp_org(rn,rn$sp.id,F)
    rn$spc=scod
    rn$mr=mr
    rn$mp=mp
    rn$nsp=nsp
    return(rn)
  }
  # return(rn)
  m=mp-mr
  
  rsp <- which(sapply(1:length(m[,1]),function(i) any(m[i,]!=0))) # especies reactivas
  csp <- which(sapply(1:length(m[,1]),function(i) all(m[i,]<=0))) # especies consumidas o no reactivas
  rcsp <- intersect(rsp,csp) # especies consumidas y reactivas
  
  if(length(rcsp!=0)){
    cm <- matrix(0,length(scod),length(rcsp))
    mr <- cbind(mr,cm)
    
    for(i in 1:length(rcsp)){
      cm[rcsp[i],i]=1
      rn$reac <- list.append(rn$reac,list(reactants=NULL,products=rcsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
    }
    mp <- cbind(mp,cm)
    
    
  }
  
  
  rn$spc=scod
  rn$mr=mr
  rn$mp=mp
  
  nsp <- rn.linp_org(rn,rn$sp.id,F)
  if(length(nsp)==0) {
    return(rn)
  }
  
  cm <- matrix(0,length(scod),length(nsp))
  rn$mr <- cbind(rn$mr,cm)
  
  for(i in length(nsp)){
    cm[nsp[i],i]=1
    rn$reac <- list.append(rn$reac,list(reactants=NULL,products=nsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
  }
  rn$mp <- cbind(rn$mp,cm)
  
  return(rn)
}

b_struct <- function(){
  
  b_levels <- sapply(crn_r$p_b,FUN = length)
  b_features <- function(i){
    r <- list()
    r$reactions.partition <- crn_r$r_p[[i]]
    r$species.partition <- crn_r$sp_p[[i]]
    r$level=length(crn_r$p_b[[i]])
    r$basic.contained <- crn_r$p_b[[i]]
    r$basic.level.contianied <- table(b_levels[crn_r$p_b[[i]]])
    r$reactions.basic <- crn_r$r_b[[i]]
    r$species.basic <- crn_r$sp_b[[i]]
    dcomp <- rn.dcom(rn,r$species.basic)
    r$species.overproduce <- dcomp$opsp
    r$species.fragilecycle <- dcomp$fsp
    r$species.cycle.correspondence <- dcomp$ec
    return(r)
  }
  return(sapply(1:length(crn_r$p_b),FUN=b_features))
}

org_struct <- function(tot_org=T){
  
  b_levels <- sapply(crn_r$p_b,FUN = length)
  b_features <- function(i){
    r <- list()
    r$reactions.partition <- crn_r$r_p[[i]]
    r$species.partition <- crn_r$sp_p[[i]]
    r$level <- length(crn_r$p_b[[i]])
    r$basic.contained <- crn_r$p_b[[i]]
    r$basic.level.contained <- table(b_levels[crn_r$p_b[[i]]])
    r$reactions.basic <- crn_r$r_b[[i]]
    r$species.basic <- crn_r$sp_b[[i]]
    dcomp <- rn.dcom(rn,r$species.basic)
    r$species.overproduce <- dcomp$opsp
    r$species.fragilecycle <- dcomp$fsp
    r$species.cycle.correspondence <- dcomp$ec
    return(r)
  }
  out <- list()
  basics <- sapply(1:length(crn_r$p_b),FUN=b_features)
  out$basics <- basics
  
  rn_run_c_sso_gen()
  crn_r<-rn_get()
  
  org_list <- list()
  if(length(crn_r$ssmo_sp)>0) org_list <- 
    crn_r$ssmo_p[which(sapply(1:length(crn_r$ssmo_sp),function(i) rn.linp_sorg_c(crn_r,i,F)))]
  if(length(org_list)==0){
    out$species.number <- length(rn$sp.name)
    out$reaction.number <- length(rn$reac)
    out$basics.number <- length(crn_r$r_b)
    return(out)
  }
  
  out$has_inflow <- any(colSums(rn$mr)==0)
  org_cont <- lapply(org_list,function(x) which(sapply(1:length(org_list), function(i) all(org_list[[i]] %in% x))))
  org_levels  <- sapply(org_list, function(x) length(which(sapply(1:length(org_list), function(i) all(org_list[[i]] %in% x)))))
  org_terminal <- sapply(org_cont,function(x) length(which(sapply(1:length(org_cont), function(i) all(x %in% org_cont[[i]]))))==1)
  
  if(out$has_inflow){
    org_elem <- sapply(org_cont,function(x) length(x)==2)
    out$inflow_basic=which(basics["level",]==1)
  }else org_elem <- sapply(org_cont,function(x) length(x)==1)
  
  org_features <- function(i){
    r <- list()
    r$reactions.partition <- unique(unlist(crn_r$r_p[unlist(org_list[i])]))
    r$species.partition <- unique(unlist(crn_r$sp_p[unlist(org_list[i])]))
    r$level=length(org_list[[i]])
    r$basic.contained <- org_list[[i]]
    r$basic.level.contained <- table(b_levels[org_list[[i]]])
    r$reactions.org <- unique(unlist(crn_r$r_b[unlist(org_list[i])]))
    r$species.org <- unique(unlist(crn_r$sp_b[unlist(org_list[i])]))
    dcomp <- rn.dcom(rn,r$species.org)
    r$species.overproduce <- dcomp$opsp
    r$species.fragilecycle <- dcomp$fsp
    r$species.cycle.correspondence <- dcomp$ec
    r$org.contained <- org_cont[[i]]
    r$org.deep <- max(org_levels[org_cont[[i]]])
    r$org.level.contained <- table(org_levels[org_cont[[i]]])
    r$terminal <- org_terminal[i]
    r$elemental <- org_elem[i]
    return(r)
  }
  out$org <- sapply(1:length(org_list),FUN=org_features)
  out$org.number <- length(org_list)
  out$org.terminal.number <- length(which(org_terminal))
  out$org.elemenal.number <- length(which(org_elem))
  out$species.number <- length(rn$sp.name)
  out$reaction.number <- length(rn$reac)
  out$basics.number <- length(crn_r$r_b)
  # if(!tot_org) out$totalorg.neededspecies <- rn$nsp
  out$totalorg.neededspecies <- rn$nsp
  out$N=lvm$N
  out$an=lvm$an
  out$cp=lvm$cp
  out$mu=lvm$mu
  out$cm=lvm$cm
  out$am=lvm$am
  out$ind=lvm$ind
  cat("especies necesarias",rn$nsp,"\n")
  return(out)
  
}

close_struct <- function(tot_org=T){
  
  b_levels <- sapply(crn_r$p_b,FUN = length)
  b_features <- function(i){
    r <- list()
    r$reactions.partition <- crn_r$r_p[[i]]
    r$species.partition <- crn_r$sp_p[[i]]
    r$level <- length(crn_r$p_b[[i]])
    r$basic.contained <- crn_r$p_b[[i]]
    r$basic.level.contained <- table(b_levels[crn_r$p_b[[i]]])
    r$reactions.basic <- crn_r$r_b[[i]]
    r$species.basic <- crn_r$sp_b[[i]]
    dcomp <- rn.dcom(rn,r$species.basic)
    r$species.overproduce <- dcomp$opsp
    r$species.fragilecycle <- dcomp$fsp
    r$species.cycle.correspondence <- dcomp$ec
    return(r)
  }
  out <- list()
  basics <- sapply(1:length(crn_r$p_b),FUN=b_features)
  out$basics <- basics
  
  rn_run_cs_gen()
  crn_r<-rn_get()
  
  close_list <- crn_r$clos_p 
  org_list <- list()
  if(length(crn_r$clos_p)>0) org_list <- 
    crn_r$clos_p[which(sapply(1:length(crn_r$clos_sp),function(i) rn.linp_csorg_c(crn_r,i,F)))]
  if(length(close_list)==0){
    out$species.number <- length(rn$sp.name)
    out$reaction.number <- length(rn$reac)
    out$basics.number <- length(crn_r$r_b)
    return(out)
  }
  
  out$has_inflow <- any(colSums(rn$mr)==0)
  
  close_cont <- lapply(close_list,function(x) which(sapply(1:length(close_list), function(i) all(close_list[[i]] %in% x))))
  close_levels  <- sapply(close_list, function(x) length(which(sapply(1:length(close_list), function(i) all(close_list[[i]] %in% x)))))
  close_terminal <- sapply(close_cont,function(x) length(which(sapply(1:length(close_cont), function(i) all(x %in% close_cont[[i]]))))==1)
  
  org_cont <- lapply(org_list,function(x) which(sapply(1:length(org_list), function(i) all(org_list[[i]] %in% x))))
  org_levels  <- sapply(org_list, function(x) length(which(sapply(1:length(org_list), function(i) all(org_list[[i]] %in% x)))))
  org_terminal <- sapply(org_cont,function(x) length(which(sapply(1:length(org_cont), function(i) all(x %in% org_cont[[i]]))))==1)
  
  if(out$has_inflow){
    close_elem <- sapply(close_cont,function(x) length(x)==2)
    org_elem <- sapply(org_cont,function(x) length(x)==2)
    out$inflow_basic=which(basics["level",]==1)
  }else{ 
    close_elem <- sapply(close_cont,function(x) length(x)==1)
    org_elem <- sapply(org_cont,function(x) length(x)==1)
  }
  
  close_features <- function(i){
    r <- list()
    r$reactions.partition <- unique(unlist(crn_r$r_p[unlist(close_list[i])]))
    r$species.partition <- unique(unlist(crn_r$sp_p[unlist(close_list[i])]))
    r$level=length(close_list[[i]])
    r$basic.contained <- close_list[[i]]
    r$basic.level.contained <- table(b_levels[close_list[[i]]])
    r$reactions.close <- unique(unlist(crn_r$r_b[unlist(close_list[i])]))
    r$species.close <- unique(unlist(crn_r$sp_b[unlist(close_list[i])]))
    dcomp <- rn.dcom(rn,r$species.close)
    r$species.overproduce <- dcomp$opsp
    r$species.fragilecycle <- dcomp$fsp
    r$species.cycle.correspondence <- dcomp$ec
    r$close.contained <- close_cont[[i]]
    r$close.deep <- max(close_levels[close_cont[[i]]])
    r$close.level.contained <- table(close_levels[close_cont[[i]]])
    r$terminal <- close_terminal[i]
    r$elemental <- close_elem[i]
    return(r)
  }
  
  out$close <- sapply(1:length(close_list),FUN=close_features)
  out$close.number <- length(close_list)
  out$close.terminal.number <- length(which(close_terminal))
  out$close.elemenal.number <- length(which(close_elem))
  
  if(length(org_list)==0){
    out$species.number <- length(rn$sp.name)
    out$reaction.number <- length(rn$reac)
    out$basics.number <- length(crn_r$r_b)
    return(out)
  }
  
  org_features <- function(i){
    r <- list()
    r$reactions.partition <- unique(unlist(crn_r$r_p[unlist(org_list[i])]))
    r$species.partition <- unique(unlist(crn_r$sp_p[unlist(org_list[i])]))
    r$level=length(org_list[[i]])
    r$basic.contained <- org_list[[i]]
    r$basic.level.contained <- table(b_levels[org_list[[i]]])
    r$reactions.org <- unique(unlist(crn_r$r_b[unlist(org_list[i])]))
    r$species.org <- unique(unlist(crn_r$sp_b[unlist(org_list[i])]))
    dcomp <- rn.dcom(rn,r$species.org)
    r$species.overproduce <- dcomp$opsp
    r$species.fragilecycle <- dcomp$fsp
    r$species.cycle.correspondence <- dcomp$ec
    r$org.contained <- org_cont[[i]]
    r$org.deep <- max(org_levels[org_cont[[i]]])
    r$org.level.contained <- table(org_levels[org_cont[[i]]])
    r$terminal <- org_terminal[i]
    r$elemental <- org_elem[i]
    return(r)
  }
  
  out$org <- sapply(1:length(org_list),FUN=org_features)
  out$org.number <- length(org_list)
  out$org.terminal.number <- length(which(org_terminal))
  out$org.elemenal.number <- length(which(org_elem))
  out$species.number <- length(rn$sp.name)
  out$reaction.number <- length(rn$reac)
  out$basics.number <- length(crn_r$r_b)
  out$N=lvm$N
  out$an=lvm$an
  out$cp=lvm$cp
  out$mu=lvm$mu
  out$cm=lvm$cm
  out$am=lvm$am
  out$ind=lvm$ind
  
  # if(!tot_org) out$totalorg.neededspecies <- rn$nsp
  
  out$totalorg.neededspecies <- rn$nsp
  cat("especies necesarias",rn$nsp,"\n")
  return(out)
  
}


