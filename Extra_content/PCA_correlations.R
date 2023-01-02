#Make a 2D/3D PCA - correlation

x<-read.table("del_PCA_liGS.csv",sep = ';')#file input name

#Set CNV_ID as the index
rownames(x) <- x$CNV_ID

#remove ids
x$CNV_ID <- NULL

x_FP = x[x$target == 1,]
x_TP = x[x$target == 0,]

x_TP$target <- NULL
x_FP$target <- NULL

x = x_FP
#col<-c("FS1a","FS1b","FS1c","FS2a", "FS2b", "FS2c","FS3a","FS3b","FS3c")#name of the colums to be ploted
#colnames(x)<-col

x <- scale(x)
#makes the correlations
#z<-cor(x)
z<-x #no corr
#PCA distances
xx<-prcomp(z) #calculates PCA

dados<-xx$rotation
percentage<-xx$sdev^2/sum(xx$sdev^2)*100
percentage# explained variance


#plot coordinates 2D PCA
cordy=c(-0.08,0.3)
cordx=c(-0.18,0.15)

#plot 2D PCA
#biplot(xx,var.axes=FALSE,cex=c(0.6),ylabs=NULL)#, xlim=cordx, ylim=cordy)

library(ggfortify)

pca_plot <- autoplot(xx,label = FALSE, label.size = 6,size = 1.8)
#pca_plot
pca_plot + xlim(cordx) + ylim(cordy)

pca_plot <- autoplot(xx,label = TRUE, label.size = 2.8,size = 0.2)
#pca_plot
pca_plot + xlim(cordx) + ylim(cordy)

#PCA 3D
#library(rgl)
#plot3d(xx$x,xlab="PC1",ylab="PC2",zlab="PC3",type="h")
#spheres3d(xx$x, radius=0.5,col=rainbow(length(xx$x[,1])))
#grid3d(side="z", at=list(z=0))
##text3d(xx$x, text=rownames(xx$x), adj=1.3)
#text3d(xx$x, adj=1.3)