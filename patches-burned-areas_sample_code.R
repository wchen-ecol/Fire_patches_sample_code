# author: Wentao Chen, CEFE, CNRS, wentao.chen@cefe.cnrs.fr
# based on a script from Manuel campagnolo, ISA/ULisboa, mlc@isa.ulisboa.pt
# Purpose: find fire patches from Burn Area products, record patches id for each BA pixel, adapted for machines of small memory (<64Gb)
# last change: Dec., 2019


############################# TODO ###################################
#### 1. Write down the results of each f_nn_stack result instead of keeping them in the memory .........ok
#### 2. Add defensive clause to each step to avoid void block errors .........ok
#### 3. Adapt for multiple years (firecci especially) .............ok with lubridate
#### 4. Transform xyd to a raster (or / eventually shapefile)
#### 5. Reduce memory use ............. en cours
#### 
#### 

###################### R packages ######################
library(raster)
library(data.table)
library(igraph)
library(RANN)
library(doParallel) # for Windows only
library(foreach)
library(SpaDES)
library(lubridate)
require(stringr)
# library(tidyr)
library(magrittr)
library(purrr)


############################ parameters ##############################

# work_dir <- "D:\\firecci_5.1/working/"
work_dir <- "~/Fdb_project/src/RUN/sample_data/" # modify here

# res_dir <- "D:\\firecci_5.1/Results/"
res_dir <- "~/Fdb_project/src/RUN/sample_data/res/"

Ncores <- 1 # number of cores to do memory-heavy jobs
Ncores_light <- 1 # number of cores to do memory-light jobs (e.g. splitting rasters)
Nthread_dt <- 5
N<- 8 #10 # number of "horizontal" blocks (to be able to process large rasters)
DELTAD <- 5 # temporal gap: if 1, then two BA pixels will belong to the same patch if they are neighbors in space and dates differ by at most 2*D
rscale <- 1

Nminpx <- 2

# convert lon/lat into sinusoidal x/y
R<-6371007
latlong2sin<-function(lon,lat) list(x= R* (pi * lon / 180) * cos (pi * lat /180), y=R * (pi * lat /180))


K <- 16  # neighbors for nn2

########################## functions #############################

counter <- function(){
  n <- 0
  f <- function(){
    n <<- n + 1
    return(n)
  }
  return(f)
}

my_add_edges <- function(graph, edges) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  edges.orig <- ecount(graph) 
  rb <- as.integer(edges[, rbind(v1, v2)])
  on.exit( .Call(igraph:::C_R_igraph_finalizer) )
  graph <- .Call(igraph:::C_R_igraph_add_edges, graph, rb - 1)
  edges.new <- ecount(graph)
  
  if (edges.new-edges.orig != 0) {
    idx <- seq(edges.orig + 1, edges.new)
  } else {
    idx <- numeric()
  }
  rm(rb)
  gc()
  graph
}


make_graph_from_dt <- function(edges, nv, directed = FALSE) {
  # assume d is a data.table
  # directed <- F
  # system.time(names <- dt[, union(unique(v1), unique(v2))])
  g <- make_empty_graph(nv, directed=directed)
  g <- my_add_edges(g, edges)
  gc()
  g
}



create.st.patches <- function(xydate=xyd,K,DIST, delta = DELTAD[1])
{
  # input: xydate=data.table with columns x, y, dtmin and dtmax: one row per pixel
  # K=number of neighbors to explore
  # DIST=max distance in meters to be a neighbor
  
  # determine neighbors
  xydate[, dtmin := dt_lbd - delta]
  xydate[, dtmax := dt_lbd + delta]
  knn<-RANN::nn2(data = cbind(xydate$lon,xydate$lat), k=min(K + 1, nrow(xydate)), searchtype = "radius", radius=DIST)
  knn$nn.dists<-NULL
  neigh<-cbind(1:nrow(knn$nn.idx),knn$nn.idx)
  rm(knn)
  
  # create adjency matrix 
  edges<-data.table::data.table(v1 = rep(neigh[,1],ncol(neigh)-1),v2 = as.vector(neigh[,-1]))
  rm(neigh)
  gc()
  edges<-edges[v2!=0]
  edges[, dtumin1 := xydate[v1, dtmin]]
  edges[, dtumax1 := xydate[v1, dtmax]]
  edges[, dtumin2 := xydate[v2, dtmin]]
  edges[, dtumax2 := xydate[v2, dtmax]]
  
  # remove edges which time span does not intersect
  edges <- edges[dtumax1 >= dtumin2 & dtumin1 <= dtumax2]
  edges[,dtumin1:=NULL]
  edges[,dtumax1:=NULL]
  edges[,dtumin2:=NULL]
  edges[,dtumax2:=NULL]
  edges <- edges[v1 != v2, ] # May lost isolated pixels
  xydate[, c("dtmin", "dtmax") := NULL]
  gc()
  
  Gall <- make_graph_from_dt(edges = edges, nv = nrow(xydate))
  rm(edges)
  gc()
  G <- igraph::simplify(Gall)
  rm(Gall)
  gc()
  # CL<-igraph::cluster_louvain(G) # extract patches (max modularity algorithm)
  CL<-igraph::clusters(G) # extract patches (maximally connected components)
  rm(G)
  return(CL$membership)
}

######################################################


isBorder <- function(k, N = N, id_dt, nc, nr){
  # return a data.table with each row specifying whether a pixel with raster id and blk number k lies at border of a nc*nr raster, and which (t, b, l, r) border 
  dborder <- function(nc, nr){
    top <- c(1:(3 * nc))
    bottom <- (nc * (nr - 3) + 1) : (nc * nr)
    left <- c(Reduce('+', rep(nc, nr - 1), init = 1, accumulate = T),
              Reduce('+', rep(nc, nr - 1), init = 2, accumulate = T),
              Reduce('+', rep(nc, nr - 1), init = 3, accumulate = T))
    right <- c(Reduce('+', rep(nc, nr - 1), init = nc - 2, accumulate = T),
               Reduce('+', rep(nc, nr - 1), init = nc - 1, accumulate = T),
               Reduce('+', rep(nc, nr - 1), init = nc, accumulate = T))
    return(data.table(id = c(top, bottom, left, right),
                      direction = c(rep("t", length(top)),
                                    rep("b", length(bottom)),
                                    rep("l", length(left)),
                                    rep("r", length(right)))))
  }
  dborderN <- function(nc, nr){
    top <- Reduce('+', rep(N, N -1), init = N, accumulate = T)
    bottom <- Reduce('+', rep(N, N -1), init = 1, accumulate = T)
    left <- 1:N
    right <- (N * N):(N*(N - 1) + 1)
    return(data.table(id = c(top, bottom, left, right),
                             direction = c(rep("t", length(top)),
                                           rep("b", length(bottom)),
                                           rep("l", length(left)),
                                           rep("r", length(right)))))
  }
  
  bds <- dborder(nc, nr)
  bds_to_retain <- bds[!(direction %chin% dborderN(N,N)[id == k, direction]), ]
  setkey(bds_to_retain, 'id')
  id_dt[, ':='(isbd_t = id %in% bds_to_retain[direction == 't', id],
               isbd_b = id %in% bds_to_retain[direction == 'b', id],
               isbd_l = id %in% bds_to_retain[direction == 'l', id],
               isbd_r = id %in% bds_to_retain[direction == 'r', id])]
  id_dt[, isbd := isbd_t | isbd_b | isbd_l | isbd_r]
  id_dt[, c("isbd_t", "isbd_b", 'isbd_l', 'isbd_r') := NULL]
  return(id_dt)
}

#########################################

f_nn_stack <- function(rpart, delta = DELTAD[1]){
  # main function: to be applied to the i-th  block 
  # it returns a data.table xyd with column "membership" recording patch ids
  #  convert raster to data.table & remove unburned and unclassified (dt>0)
  make_xyd <- function(rblk, Nlayer){
    vr <- rblk[]
    cells <- cellFromRowCol(rblk, c(1, nrow(rblk)), c(1, ncol(rblk)))
    coordsr <- xyFromCell(rblk, 1:ncell(rblk))
    xyd <- data.table(x=coordsr[,1], y=coordsr[,2], dt=vr) 
    xyd[, id := 1:ncell(rblk)]
    xyd <- xyd[dt>0] 
    colnames(xyd) <- c("lon","lat","dt", "id")
    

    system.time(xyd[, dt_lbd := as.Date(dt, origin = as.Date(paste(YEAR[Nlayer], "01-01", sep = "-")))])
    xyd[, Nlayer := Nlayer]
    return(xyd)
  }
  xyd <- rbindlist(lapply(1:length(rpart), function(i)make_xyd(rpart[[i]], i)))
  # determine patches
  if(nrow(xyd) == 0){
    xyd[, membership := NA]
    return(xyd)
  }else{
    xyd$membership <- create.st.patches(xydate=xyd,K=K,DIST=d, delta = delta) #
    return(xyd)
  }
}



###################################### main ################################


# getDTthreads()
# setDTthreads(3)
setwd(work_dir)

## in bash " for FN in `ls |grep 'h41v20'`; do tar -xf $FN -C /mnt/c/Users/WCHEN/Documents/Fdb_project/src/RUN/; done"

# fn <- list.files(".", pattern = paste(LOCATION, "\\-fv1\\.1\\-JD\\.tif$", sep = "")) # for sentinel2
# timer <- data.table(idx = 0, loc = "", Tstart = Sys.time(), Tend = Sys.time(), Tdiff = 0.)
# for(idx in c(1:length(LOCATIONs))){ 
  # tstart <- Sys.time()
  task_number <- 1
  LOCATION <- "h39v15" # LOCATIONs[idx]
  fn_pattern <- paste(LOCATION, "\\-fv1\\.1\\-JD\\.tif$", sep = "") 
  
  fn <- list.files(".", pattern = fn_pattern) # for firecci 4.1
  
  # if(length(fn) < 1){
  #   next
  # }
  
  TIME <- stringr::str_sub(fn, start = 1, end = 8) # Africa cci specific
  TIME <- ymd(TIME)
  
  YEAR <- as.integer(stringr::str_sub(TIME, 1, 4))
  MONTH <- as.integer(stringr::str_sub(TIME, 5, 6))
  # MONTH
  
  rstack <- stack(fn)
  rstack[[1]]@extent == rstack[[2]]@extent
  
  
  
  if(rscale > 1){
    rstack <- crop(rstack,extent(rstack)/rscale) # rscale^2 times fewer pixels
  }
  
  res<-res(rstack)
  ext<-as.vector(extent(rstack))
  
  deltalon<-res[1]; deltalat<-res[2]
  
  d <- sqrt(res[1]^2+res[2]^2) # d~30 meters
  if(2 * res[1] <= d | 2 * res[2] <= d)stop("irregular resolution")
  d <- d + 0.5 * (d - min(res[1], res[2]))
  d
  
  ###########
  
  task_tmpdir <- paste("TMP", task_number, "/", sep = "")
  dir.create(task_tmpdir)
  res_dir <- paste("./Res", LOCATION, "/", sep = "_")
  dir.create(res_dir)
  list_r <- list()
  
  system.time(for(i in 1:nlayers(rstack)){
    list_r[i] <- rstack[[i]]
  })
  
  cl <- makeCluster(Ncores_light)
  
  registerDoParallel(cl)
  
  
  ### split each layer in the stack into N*N blocks, write to hard disk in order to reduce memory use
  system.time(r_splitted <- foreach(lr = list_r, .packages = c('raster', 'SpaDES'))%dopar%{
      tmp <- splitRaster(lr, nx = N, ny = N)
      fns <- c()
      nc <- c()
      nr <- c()
      for(i in 1:length(tmp)){
        fn <- paste(task_tmpdir, i, tmp[[i]]@data@names, "scale1_", rscale, '.tif', sep = '_')
        writeRaster(tmp[[i]], filename = fn, format = 'GTiff', overwrite = T)
        fns <- c(fns, fn)
        nc <- c(nc, tmp[[i]]@ncols)
        nr <- c(nr, tmp[[i]]@nrows)
      }
      list(fn = fns, nc = nc, nr = nr)
  })
  
  stopCluster(cl)
  
  gc()
  
  length(r_splitted)
  
  ## regroup blocks of the same location (but with different periods) by file names
  system.time(r_part <- lapply(1:(N * N), 
                               function(x)list(x, sapply(r_splitted, 
                                                 function(y)y$fn[x]))))
  # r_part
  
  ## record the dimensions of each block, assuming blocks of the same location but of different periods have exactly the same dimensions
  ncnr <- data.table(nc = r_splitted[[1]]$nc, nr = r_splitted[[1]]$nr)
  # ncnr
  
  
  ## compute patches for each block
  
  FN_dt <- paste(task_tmpdir, YEAR[1], LOCATION, "dt", sep = '_')
  
  cl <- makeCluster(Ncores)
  registerDoParallel(cl)
  
  system.time(lista <- foreach(rp_names = r_part, .packages = c("RANN", "data.table", "raster", "igraph"))%dopar%{
    setDTthreads(3)
    rs <- stack(rp_names[[2]]) # load blocks from hard disk
    rp <- list()
    for(i in 1:nlayers(rs)){
      rp[i] <- rs[[i]]
    }
    res <- f_nn_stack(rp)
    res <- unique(res)
    setkeyv(res, c('id', 'dt_lbd'))
    res[, blk := rp_names[[1]]]
    bd <- isBorder(rp_names[[1]], N, res[, .(id, dt_lbd)], ncnr$nc[i], ncnr$nr[i])
    setkeyv(bd, c('id', 'dt_lbd'))
    res <- res[bd, ] # determine the need to regroup to avoid border disruption, merge bd and res by keys: id, dt_lbd
    FN_dt_i <- paste(FN_dt, formatC(rp_names[[1]], width = 3, flag = "0"), "nonbd", ".RData", sep = "_")
    FN_dt_i_bd <- paste(FN_dt, formatC(rp_names[[1]], width = 3, flag = "0"), "isbd", ".RData", sep = "_")
    bd_patches <- unique(res[isbd == T, membership])
    res_bd <- res[membership %in% bd_patches, ]
    res_nonbd <- res[!(membership %in% bd_patches), ]
    saveRDS(res_nonbd, file = FN_dt_i, compress = F)
    saveRDS(res_bd, file = FN_dt_i_bd, compress = F)
    list(Fn = list(FN_dt_i, FN_dt_i_bd, res_nonbd[, .N], res_bd[, .N]), 
         Nid = res_nonbd[, ifelse(.N > 0, max(id), 0)], 
         Npatches = res_nonbd[, ifelse(.N > 0, max(membership), 0)])
    })
  
  stopCluster(cl)
  gc()
  # Ncores <- 3
  Nxyd <- sum(sapply(lista, function(x)x[[1]][[3]])) # total number of pixels
  Nxyd_bd <- sum(sapply(lista, function(x)x[[1]][[4]])) # total number of pixels of border patches
  FN_nonbd <- sapply(lista, function(x)x[[1]][[1]]) # file names of pixels of non-border patches
  FN_isbd <- sapply(lista, function(x)x[[1]][[2]]) # file names of pixels of border patches
  Nid <- Reduce("+", sapply(lista, function(x)x[[2]]), init = 0, accumulate = T) # accumulated numbers of ids
  Npatches<- Reduce("+", sapply(lista, function(x)x[[3]]), init = 0, accumulate = T)
  
  
  if(Nxyd < 1){ # save non border data
    xyd_reduced <- data.table(IDmap= NA, id = NA, lon = NA, lat = NA, Nlayer = NA, blk = NA, dt_lbd = NA)
    xyd_reduced[, membership_combined := NA]
    system.time(saveRDS(xyd_reduced, file = paste(res_dir, YEAR[1], LOCATION, "BUFFER", DELTAD, "N", N, "Npx", Nminpx, "scale", rscale, "nonbd.void.RData.gzip", sep = "_"), compress = T))
  }else{
    cl <- makeCluster(Ncores)
    registerDoParallel(cl)
    
    system.time(foreach(idx = 1:length(FN_nonbd), 
                        .packages = c("data.table"))%dopar%{
                          tmp <- readRDS(FN_nonbd[idx])
                          tmp[, IDmap := id + Nid[idx], by = blk]
                          tmp[, membership_combined := membership + Npatches[blk], by = blk]
                          xyd_reduced <- tmp[, .(IDmap, id, lon, lat, Nlayer, blk, dt_lbd, membership_combined)]
                          system.time(saveRDS(xyd_reduced, file = paste(res_dir, YEAR[1], LOCATION, "BUFFER", DELTAD, "N", N, "Npx", Nminpx, "scale", rscale, "blk", formatC(idx, width = 3, flag = "0"), "nonbd.RData.gzip", sep = "_"), compress = T))
                        })
    
    stopCluster(cl)
    gc()
  }
  
  system.time(xyd_bd <- FN_isbd %>% map(readRDS) %>% rbindlist())
  
  # setkeyv(xyd_bd, c("id", "dt_lbd"))
  
  xyd_bd <- unique(xyd_bd)
  dim(xyd_bd)
  if(nrow(xyd_bd) > 0){
    system.time(xyd_bd$membership_bd <- create.st.patches(xyd_bd, K, d, DELTAD[1])) # find patches with all pixels of border patches, slowest part, consumes a large amount of memory
  }else{
    xyd_bd[, membership_bd := NA]
  }
  
  # xyd_bd
  
  Ncomp_bd <- uniqueN(xyd_bd$membership_bd)
  setkey(xyd_bd, "membership_bd")
  xyd_bd[, membership_bd := membership_bd + max(Npatches)]
 
  xyd_bd[, membership_combined := membership_bd]
  xyd_bd[, membership_bd := NULL]

  xyd_bd[, IDmap := id + Nid[blk], by = blk]

  

  uniqueN(xyd_bd$membership_combined)
  
  ## remove patches smaller than Nminpx pixels (not done for now)
  
  # xyd[, pxcount := .N, by = membership_combined] # etc.
  # 
  ###### plot frequency~size relationship
  tmp <- table(xyd_bd[, .N, by = membership_combined]$N)
  freqs <- data.table(f = as.integer(tmp), size = as.integer(names(tmp)))
  rm(tmp)
  freqs
  m <- freqs[, summary(lm(log(f) ~ log(size)))]
  freqs[, plot(log(f) ~ log(size), xlab = paste("log(size)", LOCATION, sep = " "), ylab = paste("log(f) = ", round(m$coefficients[2, 1], 3), "log(size) + ", round(m$coefficients[1, 1], 3), sep = " "))]
  rm(freqs)

  ## save results as .RData
  xyd_reduced <- xyd_bd[, .(IDmap, id, lon, lat, Nlayer, blk, dt_lbd, membership_combined)]
  # setnames(xyd_reduced, "V6", "date")
  xyd_reduced
  # system.time(fwrite(xyd, file=paste(res_dir, YEAR[1], "ESACCI-L3S_FIRE-BA-MSI-AREA", LOCATION, "N", N, "Npx", Nminpx, "scale", rscale, ".txt", sep = "_")))
  system.time(saveRDS(xyd_reduced, file = paste(res_dir, YEAR[1], LOCATION, "BUFFER", DELTAD, "N", N, "Npx", Nminpx, "scale", rscale, "isbd.RData.gzip", sep = "_"), compress = T))
  rm(xyd_reduced)
  # rm(xyd_bd)
  gc()
  
  
  
  
  ####################################################################### remove temporal files
  unlink(substr(task_tmpdir, 1, str_length(task_tmpdir) - 1), recursive = T)
  # tend <- Sys.time()
  # timer <- rbindlist(list(timer, data.table(idx = idx, loc = LOCATION, Tstart = tstart, Tend = tend, Tdiff = time_length(tend - tstart) / 3600)))
# }
# timer



########## plot border patches #############
Ncomp_c <- uniqueN(xyd_bd$membership_combined)
setkey(xyd_bd, "membership_combined")
xyd_bd

dvd <- max(log(xyd_bd[, .N, by = membership_combined] + 1))

rcol <- randomcoloR::randomColor(Ncomp_c)

ct <- counter()

# pdf(file = 'test_scale30_h41v20_lonlat_reduced.pdf', height = 25, width = 25)
plot(xyd_bd$lon, xyd_bd$lat, pch="")
for (m in unique(xyd_bd$membership_combined)){
  # for (m in 1215:1219){
  tmp <- xyd_bd[membership_combined == m, ]
  # initial <- tmp[dtmin == min(dtmin), ]
  points(tmp$lon,
         tmp$lat,
         col = rcol[ct()],
         pch=".",cex=2)
  # text(tmp[, mean(x)],
  #      tmp[, mean(y)],
  #      labels = m, cex = 0.5 * log(NROW(tmp) + 1)/dvd)
#   points(initial$x,
#          initial$y,
#          col="red",
#          pch="X",cex=0.3)
}
# dev.off()


##
# mindiff <- function(dt){
#   min(diff(sort(unique(dt))))
# }


