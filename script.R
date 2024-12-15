library(vegan)


data<-read.table("termites_larvae.txt", sep="\t", h=T)

#trasform data to make them suitable for Bray-Curtis
for(c in 4:ncol(data)){
	if(min(data[,c])<0){
		data[,c]<-data[,c]-min(data[,c])
	}
}

#Calculate Bray-Curtis distances
dist<-as.matrix(vegdist(data[,4:ncol(data)]))

#Extraxct termite workers pairs
workers<-which(data$species==1)
distw<-as.vector(as.dist(dist[workers,workers]))

#Extract colony membership
colony<-data[workers,2]

#Obtain same colony vector
same<-as.vector(dist(colony))
same[which(same >0)]<-1
wilcox.test(distw,same)
boxplot(distw,same)

#Extract larvae-termites comparisons
larvae<-which(data$species==2)
dista<-dist[larvae,workers]
memb<-c(2,1,3)
same1<-NULL
otherl<-NULL

#Obtain same colony mambership
for(c in 1:3){
	sam<-which(colony==memb[c])
	same1<-c(same1,dista[c,sam])
	otherl<-c(otherl,dista[c,(c(1:ncol(dista))[-sam])])
}

#prepare the final dataset
all<-c(distw,same1,otherl)
allfact<-c(same,rep(2,length(same1)),rep(3,length(otherl)))

boxplot(all~allfact)

#allfact 0=nestmate workers, 1=non nestmate workers, 2, larvae-nestmate workers, 3, larvae non nestmate workers

#Make the KW test with 1000 resampling
kw<-kruskal.test(all~allfact)
stat<-kw$statistic
statr<-rep(NA,1000)
for(res in 1:1000){
	allfact1<-allfact[sample(1:length(allfact))]
	kw1<-kruskal.test(all~allfact1)
	statr[res]<-kw1$statistic
}
length(which(statr>stat))/1000

#Select the pairs for pairwise comparisons

pairs<-c(1,3)
all2<-all[which(allfact %in% pairs)]
allfact2<-allfact[which(allfact %in% pairs)]

kw<-wilcox.test(all2~allfact2)
stat<-kw$p.value
statr<-rep(NA,1000)
for(res in 1:1000){
	allfact1<-allfact2[sample(1:length(allfact2))]
	kw1<-wilcox.test(all2~allfact1)
	statr[res]<-kw1$p.value
}
length(which(statr>stat))/1000


#Make the MDS plot
colours<-c("red","red","blue","blue")
pch<-c(1,16,1,16)
mds<-metaMDS(dist)
mds$points
plot(mds$points,col=colours[data[,2]],pch=pch[data[,3]])





