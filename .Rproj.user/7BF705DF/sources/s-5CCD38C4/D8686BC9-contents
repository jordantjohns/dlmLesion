#' Function for plotting data (predicted and original)
#'
#' Needs two inputs, a data frame dat_b, and a vector of vids (voxels that are going to be highlighted)
#' plot_dlm is a function for just one voxel
#' plot_dlm2 is a function for two voxels
#'
#' @param dat_b is a data frame with the following columns/variables
#' lid (lesion id)
#' vid (voxel id)
#' xind,yind,zind (x,y,z coordinates for voxel location in lesion/brain)
#' y (voxel intensity)
#' source (character vector, indicating source of voxel intensity data - with two options)
#' "Standardized" or "Predicted"
#' stime (time of visits, can be interpolated times or actual times - numeric)
#'
#' @param vids is one voxel id from the data (the voxel that will be plotted)
#'
#' @return plot object (grob) from grid.arrange, also outputs the plot to window
#' @export
#'

plot_dlm<-function(dat_b,vids){
  library(ggplot2)
  library(gridExtra)
  #checks


  ttt1<- dat_b %>% filter(vid==vids)
  lid1<-ttt1$lid[1]
  xind1<-ttt1$xind[1]
  #id1<-ttt1$id[1]

  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }

  #plot A
  less_dat<-dat_b %>% filter(lid==lid1 &xind==xind1 & (!is.na(y)) ) %>%
    arrange(vid,source,stime)

  col_norm<-"gray25"
  col_pred<-"red"

  bground_col='transparent'
  lines_p<-ggplot(less_dat,aes(x=stime,y=y,group=interaction(vid,source))) +
    geom_line(aes(color=source),alpha=.15) + scale_colour_manual("",values=c(col_pred,col_norm))
  lines_p2 <- lines_p + geom_line(data=filter(less_dat,vid==vids), aes(x=stime,y=y,group=source,colour=source))
  lines_p2 <- lines_p2+ labs(x="Time from Lesion Incidence (days)",y="Intensity") + ggtitle(paste0("A: Voxel Intensity Trajectories") ) +
    theme(legend.position = "none")

  #slice 1 - B
  yinda<-filter(less_dat,source=="Standardized" & vid == vids)$yind[1]
  zinda<-filter(less_dat,source=="Standardized" & vid == vids)$zind[1]

  ann_text <- data.frame(stime = min(less_dat$stime),yind = yinda-.75, zind=zinda-.45,source="Standardized",lid=lid1,xind=xind1,vid=vids)


  sliceX_p1=ggplot(filter(less_dat,source=="Standardized"), aes(x=yind,y=zind)) + geom_tile(aes(fill = y))
  sliceX_p1<-sliceX_p1 + facet_grid(source~stime) + scale_fill_gradient2(guide=guide_colorbar("FLAIR Intensity", title.position="top", barwidth = 10), low="midnightblue",
                                                                         mid="white", high="firebrick", midpoint = 2.5,
                                                                         limits=c(min(less_dat$y,na.rm=T),max(less_dat$y,na.rm=T)))+
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),strip.background = element_blank(),
          strip.text.y = element_blank(),legend.position="bottom",legend.direction="horizontal") +
    labs(y="Standardized", x="") +
    ggtitle(paste0("B: L", lid1," Slice, Standardized Data")) +
    annotate('segment',x=yinda-1,xend=yinda,y=zinda-.25,yend=zinda-.25, colour="black", size=.5, alpha=.9, arrow=arrow(length=unit(0.05,"cm"))) #+

  leg_p1<-g_legend(sliceX_p1)
  sliceX_p1=sliceX_p1 + theme(legend.position="none")

  less_2<-less_dat %>% filter(source=="Predicted")

  #slice 1 - C
  ann_text <- data.frame(stime = min(less_2$stime,na.rm=T),yind = yinda-.75, zind=zinda-.45,source="Predicted",lesionid=lid1,xind=xind1,vid=vids)


  sliceX_p2=ggplot(filter(less_dat,source=="Predicted"), aes(x=yind,y=zind)) + geom_tile(aes(fill = y))
  sliceX_p2<-sliceX_p2 + facet_grid(source~stime) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), strip.background = element_blank(),
          strip.text.y = element_blank(), axis.title.x=element_text(hjust=0,size=12)) +
    scale_fill_gradient2(low="midnightblue",
                         mid="white", high="firebrick", midpoint = 2.5,
                         limits=c(min(less_dat$y,na.rm=T),max(less_dat$y,na.rm=T)))+
    labs(y="Predicted", x="Time from incidence (days)", fill="FLAIR (P)") +
    ggtitle(paste0("C: L", lid1," Slice, Predicted Data")) +
    annotate('segment',x=yinda-1,xend=yinda,y=zinda-.25,yend=zinda-.25, colour="black", size=.5, alpha=.9, arrow=arrow(length=unit(0.05,"cm"))) #+

  sliceX_p2=sliceX_p2 + theme(legend.position="none")

  #good
  laymtx <- rbind(c(1,1,1,1),
                  c(1,1,1,1),
                  c(1,1,1,1),
                  c(1,1,1,1),
                  c(2,2,2,2),
                  c(2,2,2,2),
                  c(4,4,3,3),
                  c(4,4,3,3))

  slices<-list(lines_p2,sliceX_p1,sliceX_p2,leg_p1)
  plotres<-grid.arrange(grobs=slices,layout_matrix = laymtx)
  return(plotres)
}

