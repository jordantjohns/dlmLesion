#' Predicts intensities for new voxels
#'
#'
#' @param newdata_df is a data frame with new voxel information, including variables for all terms in model
#' can be one or multiple voxels, but must be a dataframe
#' @param curr_models is a list of linear model objects, comes directly from output from dlm command (dlm_model)
#' @param ptimes is vector of prediction times, use for column names for new data predictions
#' @return pred_obs-a matrix of prediction intensities at the times given in the input ptimes
#' @export
#'

dlm_pred<-function(newdata_df,curr_models,ptimes){
  #checks

  #predict new data at each prediction time (specified by ptimes)
  pred_obs<-matrix(NA,ncol=length(curr_models),nrow=nrow(newdata_df))

  for(i in 1:length(curr_models)){
    pred_obs[,i] <- as.vector(predict(curr_models[[i]], newdata = newdata_df))

  }
  colnames(pred_obs)<-paste0("t",ptimes)
  return(pred_obs)
}




