#' Interpolation functions
#' time_interpolate is a helper function for TimeInterpolateByGroup
#' that operates on each of the groups. In the input to this function,
#' the grouping_variable column of the data frame should be single-valued.
#' The function returns a (probably longer) data frame, with estimated
#' values for the times specified in the output_times array.

time_interpolate <- function(data_frame,
                             grouping_variable,
                             time_var,
                             output_times) {
  input_times <- data_frame[, time_var]
  exclude_vars <- c(time_var, grouping_variable)
  value_vars <- setdiff(colnames(data_frame), exclude_vars)
  output_df <- data.frame(rep(data_frame[1,grouping_variable], length(output_times)), output_times)
  colnames(output_df) <- c(grouping_variable, time_var)
  for (value_var in value_vars) {
    output_df[,value_var] <- approx(input_times, data_frame[, value_var], output_times)$y
  }
  return(output_df)
}

time_extrapolate <- function(data_frame,
                             grouping_variable,
                             time_var,
                             output_times) {
  input_times <- data_frame[, time_var]
  exclude_vars <- c(time_var, grouping_variable)
  value_vars <- setdiff(colnames(data_frame), exclude_vars)
  output_df <- data.frame(rep(data_frame[1,grouping_variable], length(output_times)), output_times)
  colnames(output_df) <- c(grouping_variable, time_var)
  for (value_var in value_vars) {
    output_df[,value_var] <- approxExtrap(input_times, data_frame[, value_var], output_times)$y
  }
  return(output_df)
}

tibg <- function(data_frame, grouping_variable, time_variable,output_times=c(-270,-240,-210,-180,-150,-120,-90,-60,-30,0,30,60)){
  #min_time <- min(data_frame[, time_variable])
  #max_time <- max(data_frame[, time_variable])
  #output_times <- seq(from=min_time, to=max_time, by=TimeInterval)
  #output_times <- c(-120,-90,-60,-30,0,30,60,90)
  ddply(data_frame,
        grouping_variable,
        time_interpolate,
        grouping_variable=grouping_variable,
        time_var=time_variable,
        output_times=output_times)
}


tebg <- function(data_frame, grouping_variable, time_variable,output_times=c(-270,-240,-210,-180,-150,-120,-90,-60,-30,0,30,60)){
  #min_time <- min(data_frame[, time_variable])
  #max_time <- max(data_frame[, time_variable])
  #output_times <- seq(from=min_time, to=max_time, by=TimeInterval)

  ddply(data_frame,
        grouping_variable,
        time_extrapolate,
        grouping_variable=grouping_variable,
        time_var=time_variable,
        output_times=output_times)
}


#' The main function takes data from MS lesions and interpolates the times to a set grid of times
#' @return original data (df) with the following columns/variables:
#' itime = interpolated time points
#' y_i = interpolated intensity values
#'
#' The following variables must be provided:
#' @param xdat = data_frame including variables vid, stime, y_std
#' @param times = times at which we interpolate the intensities
#' @param tstar = "current time" (divides historical times and future prediction times)
#' Interpolation is done linearly before tstar and flat interpolation after tstar
#' @export
#'


interp_times<-function(xdat, times, tstar){
  #Checks

  if(!is.data.frame(xdat)) stop("xdat must be a data frame")
  require(Hmisc)

  #historical time
  t_h<-times[which(times<=tstar)]
  t_future<-times[which(times>tstar)]
  xdat2<-xdat %>% select(vid,stime,y_std)
  t_first=min(t_h)

  mid_future<-(t_future[-length(t_future)] +t_future[-1])/2
  end1 = t_future[1]-(mid_future[1]-t_future[1])
  end2 = t_future[length(t_future)]+(t_future[length(t_future)]-mid_future[length(mid_future)])
  mid_future=c(end1,mid_future,end2)

  #interpolate times
  interp_temp<-tibg(xdat2,"vid","stime",t_h)

  #extrapolate forward if needed (e.g. interpolating tstar)
    #(if there are NAs prior to  first visit - locb)
  interp_temp2<-interp_temp %>% filter(!(is.na(y_std) & stime>0))
  interp_temp3<-tebg(interp_temp2,"vid", "stime",t_h) #interpolates, but only extrapolates forward, *think this means that values missing from pre_inc are held as constant,

  #flat interpolation of response data
  flat_interp<-function(df_r, index){
    temp_resp<-df_r %>% filter(stime >= mid_future[index] & stime < mid_future[index+1]) %>%
      group_by(vid)%>%
      dplyr::summarise(y_i=mean(y_std,na.rm=T))%>%
      mutate(itime = t_future[index]) %>%
      select(vid,itime,y_i)
    return(temp_resp)
  }
  interp_resp<-flat_interp(xdat2,1)
  for(i in 2:length(t_future)){
    temp11<-flat_interp(xdat2,i)
    interp_resp<-interp_resp %>% bind_rows(temp11) %>%
      arrange(vid,itime)
  }

  interp_dat<- interp_temp3 %>% select(vid,itime=stime,y_i=y_std) %>%
    bind_rows(interp_resp) %>%
    arrange(vid,itime)

  #delete variables created
  rm(interp_resp,interp_temp,interp_temp3,interp_temp2,xdat2,end1,end2,mid_future,t_first,t_future,t_h)
  #return data frame
  return(interp_dat)
}



