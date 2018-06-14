
#################################
## R-code for Performance Plot ## 
#################################

CVdat = read.delim("CVerrors.txt",as.is=T)

minError = min(CVdat$CVerror)
minThres = CVdat$Threshold[which(CVdat$CVerror==minError)][1]
genesurv = CVdat$EdgesSelected[which(CVdat$CVerror==minError)][1]

cid = grep("CVerror_", colnames(CVdat))
ll=gsub("CVerror_","",colnames(CVdat)[cid])
cc = 2:(length(cid)+1)
legend_lab = c("Overall")
for(i in 1:length(ll)) legend_lab = c(legend_lab, ll[i])

sd = sd(CVdat$CVerror[-nrow(CVdat)])
within = CVdat$Threshold[which(abs(CVdat$CVerror)< (minError+sd))]
xpt = max(within)
ypt = CVdat$CVerror[which(CVdat$Threshold==xpt)]
zpt = CVdat$EdgesSelected[which(CVdat$Threshold==xpt)]

pdf("CVplot_Penalty.pdf",height=5, width=7, useDingbats = F)
par(mar=c(4,4,7,2),mai=c(1,1,1,0.5))
plot(CVdat$Threshold,CVdat$CVerror,ylim=c(0,1),type="l",col=1,lwd=2.5, cex.axis=0.8,ylab="Mean misclassification error", xlab="Threshold")
mtext("Edges selected",side=3,line=3)
axis(side = 3,at=CVdat$Threshold, lab=CVdat$EdgesSelected,las=2,srt= 45,cex.axis=0.8)
for(i in 1:length(cid)) lines(CVdat$Threshold, CVdat[,cid[i]], col=cc[i],lty=1)
legend("topleft",legend_lab, col=c(1,cc),lty=1 ,lwd=c(2,rep(1,length(ll))),cex=0.8)

points(minThres,CVdat$CVerror[match(minThres,CVdat$Threshold)], pch=4,lwd=2,col=2)
points(xpt,ypt, pch=4,lwd=2,col=2)
text(minThres,CVdat$CVerror[match(minThres,CVdat$Threshold)],pos=1, offset=1,lab=paste0("Minimum Threshold = ",paste(round(minThres,3),collapse="\t"),"\nClassification Error = ",round(minError,3),"\nSelected Edges = ",paste(genesurv,collapse="\t")),cex=0.6,las=1)
if(xpt!=minThres) text(xpt, ypt,pos=1, offset=1,lab=paste0( "Threshold = ",paste(round(xpt,3),collapse="\t"),"\nClassification Error = ",round(ypt,3),"\nSelected Edges = ",paste(zpt,collapse="\t")),cex=0.6,las=1)
abline(h=(minError+sd),lty=2,col="grey40")
text(-0.15,(minError+sd)+0.02, pos =4, font =2,lab="1 SD above minimum threshold",cex=0.7,col="grey40")
dev.off()