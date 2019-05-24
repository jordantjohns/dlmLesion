#' Fit dynamic prediction model
#'
#' @param dat has three columns vid, itime, and y_i, same as data from output of interp_times
#' @param htimes is a vector of historical times included as covariates in model
#' @param ptimes is vector of prediction times in the future
#' @param nfit is the number of voxels to predict for model, when NULL runs all
#' @param tx is used for parallelization, it is the tx^th iteration of the code that is run
#' tx along with nfit determines which voxels are predicted from model (line 27)
#' @return list with two objects dlm_model (a list with DLM objects from each prediction time point)
#' #' and predictions, a matrix of predicted intensities for the prediction times
#' @export
#'


dlm<-function(dat,htimes,ptimes,nfit=NULL,tx=1){
  #Checks

  #set up data
  if(!is.null(nfit)){
    nx=nfit
  }else{
    nx=nrow(dat)
  }
  ypred<-matrix(nrow = nx, ncol = length(ptimes)+1)
  dat_wide<-dat %>% filter(itime %in% as.numeric(c(htimes,ptimes))) %>%
    spread(key=itime,value=y_i)
  pmodels<-list()

  #fit_mod for each voxel and timepoint
  for(i in 1:nx){
    nvox=(tx-1)*nx+i #for parallelization
    for (j in 1:length(ptimes)) {
      t_resp<-dat_wide[nvox,] %>% select(one_of(ptimes[j])) #selects the voxel to run the model on
      #get model coefficients
      if(j>length(pmodels)){
        data.mod1a <- dat_wide %>% select(one_of(ptimes[j]))
        data.mod1b <- dat_wide %>% select(-vid,-one_of(ptimes))
        data.mod1 <- cbind(data.mod1a,data.mod1b)

        #adjust names
        names(data.mod1)[1] <- "y.mod1"
        names(data.mod1)[(dim(data.mod1)[2]-length(htimes)+1):dim(data.mod1)[2]] <-
          c( paste0("x", names(data.mod1)[(dim(data.mod1)[2]-length(htimes)+1):dim(data.mod1)[2]]))

        #fit model
        pmodels[[j]] <- lm(y.mod1 ~ ., data = data.mod1)
      }
      # if(is.na(t_resp)){
      #   ypred[i, j+1] <- NA
      # }else{
        data.mod1a <- dat_wide[-nvox,] %>% select(one_of(ptimes[j]))
        data.mod1b <- dat_wide[-nvox,] %>% select(-vid,-one_of(ptimes))
        data.mod1 <- cbind(data.mod1a,data.mod1b)

        #adjust names
        names(data.mod1)[1] <- "y.mod1"
        names(data.mod1)[(dim(data.mod1)[2]-length(htimes)+1):dim(data.mod1)[2]] <-
          c( paste0("x", names(data.mod1)[(dim(data.mod1)[2]-length(htimes)+1):dim(data.mod1)[2]]))

        #fit model
        fit.mod1 <- lm(y.mod1 ~ ., data = data.mod1)

        #predict observation left out
        new.data.mod1 <- dat_wide[nvox,] %>% select(-vid,-one_of(ptimes))

        names(new.data.mod1)[(dim(new.data.mod1)[2]-length(htimes)+1):dim(new.data.mod1)[2]] <-
          c( paste0("x", names(new.data.mod1)[(dim(new.data.mod1)[2]-length(htimes)+1):dim(new.data.mod1)[2]]))
        ypred[i,j+1] <- predict(fit.mod1, newdata = new.data.mod1)
        if(j==1) ypred[i,1] <- dat_wide$vid[nvox]
      #}
    }
  }

  colnames(ypred) <- c("vid",paste0("t",ptimes))
  return(list(dlm_model=pmodels,predictions=ypred))
}

