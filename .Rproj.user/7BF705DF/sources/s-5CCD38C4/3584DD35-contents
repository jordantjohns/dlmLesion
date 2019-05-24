#' This function takes data from MS lesions and standardizes the intensities
#' based on the pre-incidence or baseline values.
#'
#' @param method There are three standaridization methods available: Voxel level, lesion level, and stabilized versions. method = "stabilized", "voxel", or "lesion"
#' @param xdat Data must be set up in long format with the following variables:
#' lid = lesion id (factor or numeric)
#' vid = voxel id (factor or numeric)
#' vtime = time variable, can be visit number or days from incidence (numeric), where 0 = time of lesion incidence.
#' y = MRI instensity
#' @return Output is original data (df) with the following new columns/variables:
#' y_std = standardized intensities at each time point
#' base_sd = baseline sd used to calculate standardized intensities (can be unique for each vid)
#' base_mean = baseline mean intensity for pre-incidence intensities for "current" vid
#' @export
#'



baseline_std<-function(xdat, method="stabilized"){
  library(plyr)
  library(tidyverse)
  library(reshape2)
  library(lubridate)
  #data checks
    #check method one of (stabilized, lesion, voxel)
  if(!is.data.frame(xdat)) stop("xdat must be a data.frame")
  if(!(method %in% c("stabilized","lesion", "voxel"))) stop("method not one of specified options")

  #grab names of original data
  begin_vars<-names(xdat)

  #calculate baseline mean/sd for voxels
  base_summary<-xdat %>% filter(stime<0) %>%
    group_by(lid,vid) %>%
    dplyr::summarise(base_mean=mean(y,na.rm=T), sd_v=sd(y,na.rm=T)) #,ndays_base=n())

  #average baseline voxel sd for lesion
  lesion_sd <- base_summary %>% group_by(lid) %>%
    dplyr::summarise(sd_l=mean(sd_v,na.rm=T))

  #stabilized sd
  xbase<- xdat %>% left_join(base_summary,by=c("lid","vid")) %>%
    left_join(lesion_sd,by=c("lid")) %>%
    mutate(sd_stabilize = sqrt(sd_v^2/2+sd_l^2/2))

  #standardized intensity
  xbase$sd <- ifelse(method=="stabilized",xbase$sd_stabilize,
                     ifelse(method=="voxel",xbase$sd_v,
                            ifelse(method=="lesion",xbase$sd_l,NA)))
  xbase_std<-xbase %>% mutate(y_std = (y-base_mean)/sd) %>%
    select(one_of(begin_vars),y_std,base_mean,base_sd = sd)

  #delete variables created

  return(xbase_std)

}
