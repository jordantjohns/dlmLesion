#' Function for plotting data (predicted and original)
#'
#' Needs four inputs, a data frame dat_b, and a vector of vids (voxels that are going to be highlighted)
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
#'
#' @param vids is a vector of length 2 with voxel ids from the data (the two voxels that will be plotted)
#' @param aadjust is vertical adjustment for arrow label in lesion slice plot for voxel in A (first voxel)
#' @param badjust is vertical adjustment for arrow label in lesion slice plot for voxel in B (second voxel)
#'
#' @return plot object (grob) from grid.arrange, also outputs the plot to window
#' @export
#'

plot_dlm2<-function(dat_b,vids,aadjust=-0.5,badjust=.5){
  library(ggplot2)
  library(gridExtra)
  #checks
  if(!is.vector(vids)) stop("vids is not a vector")
  if(length(vids) !=2) stop("vids is not vector with length 2")

  #check vids from same lesion
  ttt1<- dat_b %>% filter(vid==vids[1])
  ttt2<- dat_b %>% filter(vid == vids[2])
  lid1<-ttt1$lid[1]
  lid2<-ttt2$lid[1]
  xind1<-ttt1$xind[1]
  xind2<-ttt2$xind[1]
  if(xind1 != xind2) stop("vids not from same lesion slice")
  #id1<-ttt1$id[1]

  #function to grab legend from plot
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }

  #plot A (trajectories - predicted and standardized, highlighted 1st voxel)
  less_dat<-dat_b %>% filter(lid==lid1 &xind==xind1 & (!is.na(y)) ) %>%
    arrange(vid,source,stime)

  col_norm<-"gray25"
  col_pred<-"red"

  bground_col='transparent'
  lines_p<-ggplot(less_dat,aes(x=stime,y=y,group=interaction(vid,source))) +
    geom_line(aes(color=source),alpha=.15) + scale_colour_manual("",values=c(col_pred,col_norm))
  lines_p2 <- lines_p + geom_line(data=filter(less_dat,vid==vids[1]), aes(x=stime,y=y,group=source,colour=source))
  lines_p2 <- lines_p2+ labs(x="Time from Lesion Incidence (days)",y="Intensity") + ggtitle(paste0("A: Voxel Intensity Trajectories") ) +
    theme(legend.position = "none")

  #plot B (trajectories - predicted and standardized, highlighted 2nd voxel)
  lines_pB<-ggplot(less_dat,aes(x=stime,y=y,group=interaction(vid,source))) +
    geom_line(aes(color=source),alpha=.15) + scale_colour_manual("",values=c(col_pred,col_norm))
  lines_pB2 <- lines_pB + geom_line(data=filter(less_dat,vid==vids[2]), aes(x=stime,y=y,group=source,colour=source))
  lines_pB2 <- lines_pB2+ labs(x="Time from Lesion Incidence (days)",y="") + ggtitle(paste0("B: Voxel Intensity Trajectories") )+
    theme(legend.background = element_rect(fill = bground_col), legend.key = element_rect(color=bground_col,fill = bground_col),
          legend.position=c(0.88,0.93), legend.direction = "vertical")


  #slice 1 - C (lesion slice from corresponding lesion, standardized data)
  yinda<-filter(less_dat,source=="Standardized" & vid == vids[1])$yind[1]
  zinda<-filter(less_dat,source=="Standardized" & vid == vids[1])$zind[1]
  yindb<-filter(less_dat,source=="Standardized" & vid == vids[2])$yind[1]
  zindb<-filter(less_dat,source=="Standardized" & vid == vids[2])$zind[1]

  ann_text <- data.frame(stime = min(less_dat$stime),yind = yinda-.75, zind=zinda+aadjust,source="Standardized",lesionid=lid1,xind=xind1,vid=vids[1])
  ann_textb <- data.frame(stime = min(less_dat$stime),yind = yindb-.75, zind=zindb+badjust,source="Standardized",lesionid=lid1,xind=xind1,vid=vids[2])

  sliceX_p1=ggplot(filter(less_dat,source=="Standardized"), aes(x=yind,y=zind)) + geom_tile(aes(fill = y))
  sliceX_p1<-sliceX_p1 + facet_grid(source~stime) + scale_fill_gradient2(guide=guide_colorbar("FLAIR Intensity", title.position="top", barwidth = 10), low="midnightblue",
                                                                         mid="white", high="firebrick", midpoint = 2.5,
                                                                         limits=c(min(less_dat$y,na.rm=T),max(less_dat$y,na.rm=T)))+
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),strip.background = element_blank(),
          strip.text.y = element_blank(),legend.position="bottom",legend.direction="horizontal") +
    labs(y="Standardized", x="") +
    ggtitle(paste0("C: L", lid1," Slice, Standardized Data")) +
    annotate('segment',x=yinda-1,xend=yinda,y=zinda-.25,yend=zinda-.25, colour="black", size=.5, alpha=.9, arrow=arrow(length=unit(0.05,"cm"))) +
    annotate('segment',x=yindb-1,xend=yindb,y=zindb,yend=zindb, colour="black", size=.5, alpha=.9, arrow=arrow(length=unit(0.05,"cm")))  +
    geom_text(data = ann_text,label = "A",size=3)+
    geom_text(data = ann_textb,label = "B",size=3)

  leg_p1<-g_legend(sliceX_p1)
  sliceX_p1=sliceX_p1 + theme(legend.position="none")
  less_2<-less_dat %>% filter(source=="Predicted")

  #slice 1 - D (lesion slice from corresponding lesion, predicted data)
  ann_text <- data.frame(stime = min(less_2$stime,na.rm=T),yind = yinda-.75, zind=zinda+aadjust,source="Predicted",lid=lid1,xind=xind1,xid=vids[1])
  ann_textb <- data.frame(stime =min(less_2$stime,na.rm=T),yind = yindb-.75, zind=zindb+badjust,source="Predicted",lid=lid1,xind=xind1,xid=vids[2])

  sliceX_p2=ggplot(filter(less_dat,source=="Predicted"), aes(x=yind,y=zind)) + geom_tile(aes(fill = y))
  sliceX_p2<-sliceX_p2 + facet_grid(source~stime) +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), strip.background = element_blank(),
          strip.text.y = element_blank(), axis.title.x=element_text(hjust=0,size=12)) +
    scale_fill_gradient2(low="midnightblue",
                         mid="white", high="firebrick", midpoint = 2.5,
                         limits=c(min(less_dat$y,na.rm=T),max(less_dat$y,na.rm=T)))+
    labs(y="Predicted", x="Time from incidence (days)", fill="FLAIR (P)") +
    ggtitle(paste0("D: L", lid1," Slice, Predicted Data")) +
    annotate('segment',x=yinda-1,xend=yinda,y=zinda-.25,yend=zinda-.25, colour="black", size=.5, alpha=.9, arrow=arrow(length=unit(0.05,"cm"))) +
    annotate('segment',x=yindb-1,xend=yindb,y=zindb,yend=zindb, colour="black", size=.5, alpha=.9, arrow=arrow(length=unit(0.05,"cm"))) +
    geom_text(data = ann_text,label = "A",size=3)+
    geom_text(data = ann_textb,label = "B",size=3)
  sliceX_p2=sliceX_p2 + theme(legend.position="none")

  #Layout matrix for figures
  laymtx <- rbind(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),
                  c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),
                  c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),
                  c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),
                  c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3),
                  c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3),
                  c(5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4),
                  c(5,5,5,5,5,5,5,5,5,4,4,4,4,4,4,4,4,4))

  #arrange plots
  slices<-list(lines_p2,lines_pB2,sliceX_p1,sliceX_p2,leg_p1)
  plotres<-grid.arrange(grobs=slices,layout_matrix = laymtx)
  return(plotres)
}
