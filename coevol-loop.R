library(phytools)
library(ape)
library(phangorn)
library(nlme)
library(lme4)
library(phylolm)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggtree)
library(ggtreeExtra)
library(pheatmap)
library(igraph)
library(ggraph)
library(Oncotree)
library(ggVennDiagram)


# this function converts a species name string from the Newick format which Common Taxonomy Tree gives us into a simpler lower-case, no quotes version comparable to Kostas' dataset
convname = function(str) {
  return(tolower(gsub("\'", "", gsub("_", " ", str))))
}

# read phylogeny previously downloaded from Common Taxonomy Tree. this has previously been cleaned:
# (i) all one line; (ii) extra set of brackets wrap the whole entity; (iii) nodes names don't contain special characters
tree = read.newick("tree-for-traits-clean-mt.phy")
#my.correlation = corBrownian(phy=tree)
tree$tip.label = convname(tree$tip.label)
tree.labels = c(tree$tip.label, tree$node.label)
root = which(tree.labels=="Eukaryota")
clade.refs = Children(tree, root)
clade.names = tree.labels[clade.refs]

# read Kostas' dataset
df = read.table("MTFull22.txt", sep="\t", header=T, stringsAsFactors = TRUE)

# manually fix bugs
df$plant.growth.form[which(df$Scientific.Name=="triodia sylvina")] = NA
df$life.cycle.habit[grep("lucilia", df$Scientific.Name)] = NA

# the dataset is read as factors by default -- we'll awkwardly convert between data types
df$Scientific.Name = as.character(df$Scientific.Name)
df$Ancestry = as.character(df$Ancestry)

df$clade = ""
to.drop = c()
for(i in 1:nrow(df)) {
  tip.ref = which(tree.labels == df$Scientific.Name[i])
  if(length(tip.ref) > 0) {
    ancestors = Ancestors(tree, tip.ref)
    ancestor = ancestors[which(ancestors %in% clade.refs)]
    df$clade[i] = tree.labels[ancestor]
  } else {
    to.drop = c(to.drop, i)
  }
}

df = df[-to.drop,]

mt.df = df

both.df = mt.df[which(!is.na(mt.df$countsMT) & !is.na(mt.df$countsPT)),]
both.df$countsMTnorm = 0
both.df$countsPTnorm = 0
for(clade in unique(both.df$clade)) {
  meanMT = mean(both.df$countsMT[both.df$clade == clade])
  meanPT = mean(both.df$countsPT[both.df$clade == clade])
  both.df$countsMTnorm[both.df$clade == clade] = both.df$countsMT[both.df$clade == clade]-meanMT
  both.df$countsPTnorm[both.df$clade == clade] = both.df$countsPT[both.df$clade == clade]-meanPT
}

both.labels = which(tree$tip.label %in% both.df$Scientific.Name)
both.tree = keep.tip(tree, both.labels)

mydf2 = data.frame(label=both.df$Scientific.Name, 
                   x=both.df$countsMT, 
                   y=both.df$countsPT, 
                   xnorm = both.df$countsMTnorm, 
                   ynorm = both.df$countsPTnorm, 
                   clade=both.df$clade)
rownames(mydf2) = mydf2$label

##### loop through different subsetting experiments
#for(expt in c(1,2,3)) {
for(expt in c(1, 2, 3, 4)) {
  if(expt == 1) { 
    mydf3 = mydf2; outstr = "coevol-all.png" 
  } else if(expt == 2) { 
    to.prune = which(mydf2$label %in% c("laminaria digitata", "nitzschia alba", "prototheca bovis", "prototheca ciferrii", "epirixanthes elongata", "nepenthes ventricosa x nepenthes alata", "choreocolax polysiphoniae"))
    mydf3 = mydf2[-to.prune,] 
    mydf3 = mydf3[mydf3$clade == "Viridiplantae",]; outstr = "coevol-plants.png" 
  } else if(expt == 3) { 
    to.prune = which(mydf2$label %in% c("laminaria digitata", "nitzschia alba", "prototheca bovis", "prototheca ciferrii", "epirixanthes elongata", "nepenthes ventricosa x nepenthes alata", "choreocolax polysiphoniae"))
    mydf3 = mydf2[-to.prune,] 
    mydf3 = mydf3[mydf3$clade != "Viridiplantae",]; outstr = "coevol-not-plants.png" 
  } else if(expt == 4) { 
    to.prune = which(mydf2$label %in% c("laminaria digitata", "nitzschia alba", "prototheca bovis", "prototheca ciferrii", "epirixanthes elongata", "nepenthes ventricosa x nepenthes alata", "choreocolax polysiphoniae"))
    mydf3 = mydf2[-to.prune,]  
    outstr = "coevol-pruned.png" 
  }
  
  colnames(mydf3) = c("label", "mt", "pt", "mtnorm", "ptnorm", "clade")
  
  ##################### 
  # linear model, no clade or correlation
  mylm = lm(pt ~ mt, data=mydf3)
  g.dumb.lm = ggplot(mydf3, aes(x=mt, y=pt)) + #geom_smooth(method="lm", alpha=0.2) +
    geom_point(aes(x=mt,y=pt,color=clade)) + xlab("MT gene count") + ylab("PT gene count") +
    labs(color="Clade") +
    theme_classic()
  g.dumb.lm
  
  #####################
  # linear model with clade correction, no correlation
  mylmn = lm(ptnorm ~ mtnorm, data=mydf3)
  g.normalised.lm = ggplot(mydf3, aes(x=mtnorm, y=ptnorm)) + geom_smooth(method="lm", alpha=0.2)+ 
    geom_point(aes(x=mtnorm,y=ptnorm,color=clade)) + xlab("Δ MT gene count") + ylab("Δ PT gene count") +
    theme_classic()
  g.normalised.lm
  
  # mixed model only meaningful for all clades
  if(expt == 1 | expt == 4) {
    ctrl <- lmeControl(opt='optim');
    #####################
    # linear mixed model with random effects for clade without clade correction
    mod.lmm = lmer(pt ~ mt + (mt | clade), data = mydf3)
    mod.lmm1 = lmer(pt ~ mt + (1 | clade), data = mydf3)
    mod.a.lmm = lme(pt ~ mt, random = ~ mt | clade, control=ctrl, data = mydf3)
    mod.a.lmm1 = lme(pt ~ mt, random = ~ 1 | clade, control=ctrl, data = mydf3)
    if(AIC(mod.lmm1) < AIC(mod.lmm)) { mod.lmm = mod.lmm1; mod.a.lmm = mod.a.lmm1 }
    
    # helper functions for mixed model plots
    fit.frame = data.frame(mt=NULL, clade=NULL)
    for(clade in unique(mydf3$clade)) {
      tmp.frame = data.frame(mt = seq(from=0, to=70, by=4), clade = clade)
      fit.frame = rbind(fit.frame, tmp.frame)
    }
    fitfn = function(model) {
      return(predict(model, fit.frame, type="response"))
    }
    cifn = function(x) {
      tx = x[order(x)]
      lci = round(length(x)*0.05)
      uci = round(length(x)*0.95)
      return(c(tx[lci],tx[uci]))
    }
    
    boot<-bootMer(mod.lmm, FUN=fitfn, nsim=100, use.u=TRUE) 
    boot.se<-apply(boot$t, 2, cifn)
    fit.frame$mod.lmm.mean = predict(mod.lmm, fit.frame, type="response")
    fit.frame$mod.lmm.boot.mean = apply(boot$t, 2, mean)
    fit.frame$mod.lmm.lo = boot.se[1,]
    fit.frame$mod.lmm.hi = boot.se[2,]
    
    g.lmm = ggplot() +
      geom_ribbon(data = fit.frame, aes(x = mt, ymin = mod.lmm.lo, ymax = mod.lmm.hi, fill = factor(clade)), alpha=0.1) +
      geom_line(data = fit.frame, aes(x = mt, y = mod.lmm.boot.mean, color = factor(clade))) +
      geom_point(data = mydf3, aes(x = mt, y = pt, color = factor(clade))) +
      scale_color_brewer(palette="Spectral") + scale_fill_brewer(palette="Spectral") +
      theme_light() + xlab("MT gene count") + ylab("PT gene count") + facet_wrap(~clade) +
      theme(legend.position="none")
    
    #####################
    # mixed model with random effects for clade on clade-corrected data
    mod.nlmm = lmer(ptnorm ~ mtnorm + (mtnorm | clade), data = mydf3)
    mod.nlmm1 = lmer(ptnorm ~ mtnorm + (1 | clade), data = mydf3)
    mod.a.nlmm = lme(ptnorm ~ mtnorm, random = ~ mtnorm | clade, control=ctrl, data = mydf3)
    mod.a.nlmm1 = lme(ptnorm ~ mtnorm, random = ~ 1 | clade, control=ctrl, data = mydf3)
    if(AIC(mod.nlmm1) < AIC(mod.nlmm)) { mod.nlmm = mod.nlmm1; mod.a.nlmm = mod.a.nlmm1 }
    
    # helper functions for mixed model plots
    fit.frame = data.frame(mtnorm=NULL, clade=NULL)
    for(clade in unique(mydf3$clade)) {
      tmp.frame = data.frame(mtnorm = seq(from=-30, to=30, by=4), clade = clade)
      fit.frame = rbind(fit.frame, tmp.frame)
    }
    fitfn = function(model) {
      return(predict(model, fit.frame, type="response"))
    }
    cifn = function(x) {
      tx = x[order(x)]
      lci = round(length(x)*0.05)
      uci = round(length(x)*0.95)
      return(c(tx[lci],tx[uci]))
    }
    
    boot<-bootMer(mod.nlmm, FUN=fitfn, nsim=100) 
    boot.se<-apply(boot$t, 2, cifn)
    fit.frame$mod.nlmm.mean = predict(mod.nlmm, fit.frame, type="response")
    fit.frame$mod.nlmm.lo = boot.se[1,]
    fit.frame$mod.nlmm.hi = boot.se[2,]
    
    g.norm.lmm = ggplot() +
      geom_ribbon(data = fit.frame, aes(x = mtnorm, ymin = mod.nlmm.lo, ymax = mod.nlmm.hi, fill = factor(clade)), alpha=0.1) +
      geom_line(data = fit.frame, aes(x = mtnorm, y = mod.nlmm.mean, color = factor(clade))) +
      geom_point(data = mydf3, aes(x = mtnorm, y = ptnorm, color = factor(clade))) +
      scale_color_brewer(palette="Spectral") + scale_fill_brewer(palette="Spectral") +
      theme_light() + xlab("Δ MT gene count") + ylab("Δ PT gene count") + facet_wrap(~clade) +
      theme(legend.position="none")
  }
  
  # pagel's lambda for various observations
  phylosig(both.tree, data.matrix(mydf3)[,2], method="lambda")
  phylosig(both.tree, data.matrix(mydf3)[,3], method="lambda")
  phylosig(both.tree, data.matrix(mydf3)[,4], method="lambda")
  phylosig(both.tree, data.matrix(mydf3)[,5], method="lambda")
  
  # scatter plots using opacity to reflect phylogenetic weighting
  mydf3$phyloweight = 0
  ntips = length(both.tree$tip.label)
  plm.boot = phylolm(pt~mt, mydf3, both.tree, model="BM", boot=1000)
  plm.df = data.frame()
  for(new.x in 0:50) {
    new.pred = plm.boot$bootstrap[,1]+plm.boot$bootstrap[,2]*new.x
    plm.df = rbind(plm.df, data.frame(x=new.x, hi=quantile(new.pred, 0.975), mean=mean(new.pred), fit=plm.boot$coefficients[1]+plm.boot$coefficients[2]*new.x, lo=quantile(new.pred, 0.025)))
  }
  
  vcv = vcvPhylo(both.tree)[1:ntips,1:ntips]
  vcv.i = solve(vcv)
  for(i in 1:length(colnames(vcv))) {
    this.sum = colSums(vcv)[i]
    this.name = colnames(vcv)[i]
    mydf3$phyloweight[mydf3$label==this.name] = 1/(1+this.sum)
  }
  g.plm.opacity.1 = ggplot() +
    geom_ribbon(data=plm.df, aes(x=x,ymin=lo,ymax=hi),alpha=0.2,fill="#0000FF") +
    geom_line(data=plm.df, aes(x=x,y=mean),alpha=1,color=1) +
    geom_point(data=mydf3, aes(x=mt, y=pt, alpha=phyloweight, color=clade)) + 
    geom_text_repel(data=mydf3, aes(x=mt,y=pt,label=label,color=clade),size=2,alpha=1) +
    xlab("MT gene count") + ylab("PT gene count") + theme_classic()
  
  plmnorm.boot = phylolm(ptnorm~mtnorm, mydf3, both.tree, model="BM", boot=1000)
  plmnorm.df = data.frame()
  for(new.x in -25:25) {
    new.pred = plmnorm.boot$bootstrap[,1]+plmnorm.boot$bootstrap[,2]*new.x
    plmnorm.df = rbind(plmnorm.df, data.frame(x=new.x, hi=quantile(new.pred, 0.975), mean=mean(new.pred), fit=plmnorm.boot$coefficients[1]+plmnorm.boot$coefficients[2]*new.x, lo=quantile(new.pred, 0.025)))
  }
  
  g.plm.opacity.2 = ggplot() +
    geom_ribbon(data=plmnorm.df, aes(x=x,ymin=lo,ymax=hi),alpha=0.2,fill="#0000FF") +
    geom_line(data=plmnorm.df, aes(x=x,y=mean),alpha=1,color=1) +
    geom_point(data=mydf3, aes(x=mtnorm, y=ptnorm, alpha=phyloweight, color=clade)) + 
   # geom_text_repel(data=mydf3, aes(x=mtnorm,y=ptnorm,label=label,color=clade),size=2,alpha=1) +
    xlab("Δ MT gene count") + ylab("Δ PT gene count") + labs(color="Clade", alpha="Relative weight") +
      theme_classic()
  
  #####################
  # phylo linear model with no clade correction
  plm = phylolm(pt~mt, mydf3, both.tree, model="BM")
  
  g.tree = ggtree(both.tree, layout='circular', branch.length='none', alpha=1) %<+% mydf3 +
    geom_fruit(geom=geom_point, pwidth=0.1,
               mapping = aes(x=pt, y=label), 
               size=0.2, color="#FF0000") +
    geom_fruit(geom=geom_point, pwidth=0.1, offset=-0.09,
               mapping = aes(x=mt, y=label), 
               size=0.2, color="#0000FF") +
    theme(legend.position="none")
  
  tree.labs = c(both.tree$tip.label, both.tree$node.label)
  
  clades.interest = c("Rhodophyta", "Apicomplexa", "Streptophyta", "Chlorophyceae")
  
  # experiment with other visualisations
  ggtree(both.tree, layout='circular', branch.length='none', alpha=1) %<+% mydf3 +
    geom_fruit(geom=geom_col, pwidth=1, alpha=0.2,
               mapping = aes(x=pt, y=label), 
               size=0.2, color="#FF0000") +
    geom_fruit(geom=geom_col, pwidth=1, offset=-1, alpha=0.2,
               mapping = aes(x=mt, y=label), 
               size=0.2, color="#0000FF") +
    #geom_cladelab(node=rhodo.ref, label="Rhodophyta") +
    geom_hilight(mapping=aes(subset = node %in% which(tree.labs %in% clades.interest), fill=nodelab(both.tree, node)),
                 type = "gradient", gradient.direction = 'rt',
                 alpha = .2) #+
  #  theme(legend.position="none")
  
 
  
  ggtree(both.tree, layout='circular', branch.length='none', alpha=1) %<+% mydf3 +
    geom_fruit(geom=geom_col, pwidth=1, alpha=0.2, offset=1,
               mapping = aes(x=ptnorm, y=label), 
               size=0.2, color="#FF0000") +
    geom_fruit(geom=geom_col, pwidth=1, offset=0, alpha=0.2,
               mapping = aes(x=mtnorm, y=label), 
               size=0.2, color="#0000FF") +
    theme(legend.position="none")
  
  #####################
  # phylo linear model with clade correction
  plm.norm = phylolm(ptnorm~mtnorm, mydf3, both.tree, model="BM")
  
  g.tree.norm = ggtree(both.tree, layout='circular', branch.length='none', alpha=1) %<+% mydf3 +
    geom_fruit(geom=geom_point, pwidth=0.1,
               mapping = aes(x=ptnorm, y=label), 
               size=0.2, color="#FF0000") +
    geom_fruit(geom=geom_point, pwidth=0.1, offset=-0.09,
               mapping = aes(x=mtnorm, y=label), 
               size=0.2, color="#0000FF") +
    theme(legend.position="none")
  
  summary(mylm)
  summary(mylmn)
  summary(mod.lmm)
  summary(mod.nlmm)
  summary(plm)
  summary(plm.norm)
  
  sum.df = data.frame()
  sum.df = rbind(sum.df, data.frame(method="naive LM slope=", coeff=summary(mylm)$coefficients[2,1], p=summary(mylm)$coefficients[2,4]))
  sum.df = rbind(sum.df, data.frame(method="shifted LM slope=", coeff=summary(mylmn)$coefficients[2,1], p=summary(mylmn)$coefficients[2,4]))
  if(expt != 2) {
    sum.df = rbind(sum.df, data.frame(method="LMM slope=", coeff=summary(mod.a.lmm)$tTable[2,1], p=summary(mod.a.lmm)$tTable[2,5]))
    sum.df = rbind(sum.df, data.frame(method="shifted LMM slope=", coeff=summary(mod.a.nlmm)$tTable[2,1], p=summary(mod.a.nlmm)$tTable[2,5]))
  }
  sum.df = rbind(sum.df, data.frame(method="phylo LM slope=", coeff=summary(plm)$coefficients[2,1], p=summary(plm)$coefficients[2,4]))
  sum.df = rbind(sum.df, data.frame(method="phylo shifted LM slope=", coeff=summary(plm.norm)$coefficients[2,1], p=summary(plm.norm)$coefficients[2,4]))
   
  tstr = function(arg) {
    return(paste(c(arg[1], signif(arg[2], digits=2), " p=", signif(arg[3], digits=2)), collapse=""))
  }
  
  sf = 2
  if(expt != 2) {
    png(outstr, width=1200*sf, height=2000*sf, res=72*sf)
    ggarrange( g.dumb.lm + ggtitle(tstr(sum.df[1,])),
                  g.normalised.lm + ggtitle(tstr(sum.df[2,])),
                  g.lmm + ggtitle(tstr(sum.df[3,])),
                  g.norm.lmm + ggtitle(tstr(sum.df[4,])),
                  g.tree + ggtitle(tstr(sum.df[5,])),
                  g.tree.norm + ggtitle(tstr(sum.df[6,])),
                  g.plm.opacity.1 + ggtitle(tstr(sum.df[5,])),
                  g.plm.opacity.2 + ggtitle(tstr(sum.df[6,])),nrow=4)
  } else {
    png(outstr, width=1200*sf, height=1200*sf, res=72*sf)
    ggarrange( g.dumb.lm + ggtitle(tstr(sum.df[1,])),
                  g.normalised.lm + ggtitle(tstr(sum.df[2,])),
                  g.tree + ggtitle(tstr(sum.df[3,])),
                  g.tree.norm + ggtitle(tstr(sum.df[4,])) , 
                  g.plm.opacity.1 + ggtitle(tstr(sum.df[5,])),
                  g.plm.opacity.2 + ggtitle(tstr(sum.df[6,])),nrow=3)
  }
  dev.off()
  
  if(expt == 1) {
    g.fig.1a = g.dumb.lm + geom_text_repel(aes(label=label, color=clade), size=2, max.overlaps=10, alpha=1) #+ ggtitle(tstr(sum.df[1,]))
  }
  if(expt == 4) {
    g.fig.1b = g.plm.opacity.2 + ggtitle(tstr(sum.df[nrow(sum.df),]))
    g.fig.s1 = g.lmm + ggtitle(tstr(sum.df[3,]))
  }
}

sf = 2
png("fig-1.png", width=800*sf, height=300*sf, res=72*sf)
ggarrange( g.fig.1a, g.fig.1b+guides(colour="none"), nrow=1, labels=c("A", "B"), font.label=list(size=18))
dev.off()

png("fig-s1.png", width=400*sf, height=400*sf, res=72*sf)
ggarrange( g.fig.s1, nrow=1)
dev.off()

write.table(mydf2, "species-list.csv",  row.names=FALSE, col.names = FALSE, sep=",")

############
# look at different ecological traits

df = read.csv("species-list-info-pruned.csv")
df$broad.cellularity = "uni"
df$broad.cellularity[grepl("multi", df$cellularity, fixed=TRUE)] = "multi"
df$broad.timing = "perennial"
df$broad.timing[df$timing != "perennial"] = "shorter"
df$broad.timing[is.na(df$timing)] = NA
df$green = "no"
df$green[which(grepl("reen", df$notes) | df$herbaceous != "undefined")] = "yes"

g.algae = ggplot(df, aes(x=mtnorm,y=ptnorm,color=alga)) + 
  geom_point() + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 
g.green = ggplot(df, aes(x=mtnorm,y=ptnorm,color=green)) + 
  geom_point() + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 

g.parasite = ggplot(df, aes(x=mtnorm,y=ptnorm,color=parasite)) + 
  geom_point() + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 
g.cellularity = ggplot(df, aes(x=mtnorm,y=ptnorm,color=broad.cellularity)) + 
  geom_point() + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 
g.herbaceous = ggplot(df, aes(x=mtnorm,y=ptnorm,color=herbaceous)) + 
  geom_point() + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 
g.timing = ggplot(df, aes(x=mtnorm,y=ptnorm,color=broad.timing)) + 
  geom_point() + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 

df$Alga = df$alga
g.algae.lm = ggplot(df, aes(x=mtnorm,y=ptnorm,color=Alga,fill=Alga)) +
 geom_smooth(data=df[df$Alga=="yes",], method="lm", alpha=0.1, linewidth=0.2) + 
  geom_smooth(data=df[df$Alga=="no",], method="lm", alpha=0.1, linewidth=0.2) +
  geom_point(data = df) + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 

df$Herbaceous = df$herbaceous
df$Herbaceous[df$Herbaceous=="mixed"] = "yes"
g.herbaceous.lm = ggplot(df, aes(x=mtnorm,y=ptnorm,color=Herbaceous,fill=Herbaceous)) + 
  geom_smooth(data=df[df$Herbaceous=="yes",], method="lm", alpha=0.1, linewidth=0.2) + 
  geom_smooth(data=df[df$Herbaceous=="no",], method="lm", alpha=0.1, linewidth=0.2) +
  geom_point() + theme_classic() + xlab("Δ MT gene count") + ylab("Δ PT gene count") 

x = list(
  Unicellular = which(df$broad.cellularity=="uni"),
  Green = which(grepl("reen", df$notes) | df$herbaceous != "undefined"),
  Algae = which(df$alga=="yes"),
  #  Multicellular = which(df$broad.cellularity=="multi"),
  Herbaceous = which(df$herbaceous=="yes")
)
g.venn.1 = ggVennDiagram(x) + 
  scale_fill_gradient(low="white",high = "blue") +
  scale_color_discrete("#000000") 

x = list(
  Perennial = which(df$broad.timing=="perennial"),
  Green = which(grepl("reen", df$notes) | df$herbaceous != "undefined"),
  Algae = which(df$alga=="yes"),
  #  Multicellular = which(df$broad.cellularity=="multi"),
  Herbaceous = which(df$herbaceous=="yes")
)
g.venn.2 = ggVennDiagram(x) + 
  scale_fill_gradient(low="white",high = "blue") +
  scale_color_discrete("#000000") 
png("venn-diags.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange(g.venn.1, g.venn.2, nrow=1, labels=c("A", "B"), font.label=list(size=18))
dev.off()

#png("fig-2.png", width=400*sf, height=700*sf, res=72*sf)
#ggarrange(g.venn.1 + theme(legend.position = "none"), g.algae, g.herbaceous, nrow=3, labels=c("A", "B", "C"), font.label=list(size=18))
#dev.off()

png("fig-2.png", width=700*sf, height=300*sf, res=72*sf)
ggarrange(g.venn.1 + theme(legend.position = "none") + scale_x_continuous(expand = c(0.1, 0.1)), g.algae.lm, nrow=1, labels=c("A", "B"), font.label=list(size=18))
dev.off()

### LMMs for eco traits

ctrl <- lmeControl(opt='optim');

# unicellularity
mod.lm = lm(ptnorm ~ mtnorm, data = df)
#mod.nlmm = lmer(ptnorm ~ mtnorm + (mtnorm | broad.cellularity), data = df, method="ML")
#mod.nlmm1 = lmer(ptnorm ~ mtnorm + (1 | broad.cellularity), data = df, method="ML")
mod.a.nlmm = lme(ptnorm ~ mtnorm, random = ~ mtnorm | broad.cellularity, control=ctrl, data = df, method="ML")
mod.a.nlmm1 = lme(ptnorm ~ mtnorm, random = ~ 1 | broad.cellularity, control=ctrl, data = df, method="ML")
AIC(mod.lm, mod.a.nlmm, mod.a.nlmm1)

# green
mod.lm = lm(ptnorm ~ mtnorm, data = df)
#mod.nlmm = lmer(ptnorm ~ mtnorm + (mtnorm | broad.cellularity), data = df, method="ML")
#mod.nlmm1 = lmer(ptnorm ~ mtnorm + (1 | broad.cellularity), data = df, method="ML")
mod.a.nlmm = lme(ptnorm ~ mtnorm, random = ~ mtnorm | green, control=ctrl, data = df, method="ML")
mod.a.nlmm1 = lme(ptnorm ~ mtnorm, random = ~ 1 | green, control=ctrl, data = df, method="ML")
AIC(mod.lm, mod.a.nlmm, mod.a.nlmm1)

# algae
mod.lm = lm(ptnorm ~ mtnorm, data = df)
#mod.nlmm = lmer(ptnorm ~ mtnorm + (mtnorm | broad.cellularity), data = df, method="ML")
#mod.nlmm1 = lmer(ptnorm ~ mtnorm + (1 | broad.cellularity), data = df, method="ML")
mod.a.nlmm = lme(ptnorm ~ mtnorm, random = ~ mtnorm | alga, control=ctrl, data = df, method="ML")
mod.a.nlmm1 = lme(ptnorm ~ mtnorm, random = ~ 1 | alga, control=ctrl, data = df, method="ML")
AIC(mod.lm, mod.a.nlmm, mod.a.nlmm1)

# herbaceous
#df$herbaceous[is.na(df$herbaceous)] = "undefined"
h.df = df[df$Herbaceous == "yes" | df$Herbaceous == "no",]
mod.lm = lm(ptnorm ~ mtnorm, data = h.df)
#mod.nlmm = lmer(ptnorm ~ mtnorm + (mtnorm | herbaceous), data = df, method="ML")
#mod.nlmm1 = lmer(ptnorm ~ mtnorm + (1 | herbaceous), data = df, method="ML")
mod.a.nlmm = lme(ptnorm ~ mtnorm, random = ~ mtnorm | Herbaceous, control=ctrl, data = h.df, method="ML")
mod.a.nlmm1 = lme(ptnorm ~ mtnorm, random = ~ 1 | Herbaceous, control=ctrl, data = h.df, method="ML")
AIC(mod.lm, mod.a.nlmm, mod.a.nlmm1)

# timing
df$broad.timing[is.na(df$broad.timing)] = "undefined"
mod.lm = lm(ptnorm ~ mtnorm, data = df)
#mod.nlmm = lmer(ptnorm ~ mtnorm + (mtnorm | herbaceous), data = df, method="ML")
#mod.nlmm1 = lmer(ptnorm ~ mtnorm + (1 | herbaceous), data = df, method="ML")
mod.a.nlmm = lme(ptnorm ~ mtnorm, random = ~ mtnorm | broad.timing, control=ctrl, data = df, method="ML")
mod.a.nlmm1 = lme(ptnorm ~ mtnorm, random = ~ 1 | broad.timing, control=ctrl, data = df, method="ML")
AIC(mod.lm, mod.a.nlmm, mod.a.nlmm1)

df$green = "no"
df$green[which(grepl("reen", df$notes) | df$herbaceous != "undefined")] = "yes"
mod.lm = lm(ptnorm ~ mtnorm, data = df)
#mod.nlmm = lmer(ptnorm ~ mtnorm + (mtnorm | herbaceous), data = df, method="ML")
#mod.nlmm1 = lmer(ptnorm ~ mtnorm + (1 | herbaceous), data = df, method="ML")
mod.a.nlmm = lme(ptnorm ~ mtnorm, random = ~ mtnorm | green, control=ctrl, data = df, method="ML")
mod.a.nlmm1 = lme(ptnorm ~ mtnorm, random = ~ 1 | green, control=ctrl, data = df, method="ML")
AIC(mod.lm, mod.a.nlmm, mod.a.nlmm1)


# try and build predictive models of gene counts
#pred.lm = lm(pt ~ mt*alga + mt*broad.cellularity + mt*parasite, data=df)
#pred.lm = lm(pt ~ clade + mt*alga + mt*broad.cellularity + mt*parasite, data=df)
#pred.lm = lm(pt ~ clade*mt, data=df)
#pred.lm = lm(pt ~ clade + mt, data=df)
#pred.lm = lm(ptnorm ~ mtnorm, data=df)
#pred.rf = randomForest(pt ~ clade*mt, data=df)
#pred.tree = tree(pt ~ clade*mt, data=df)

########### explore clustering of genes across different organelles

mt.df = read.csv("mt-barcodes-manual.csv")
colnames(mt.df) = paste("MT", colnames(mt.df), sep="-")
pt.df = read.csv("pt-barcodes-manual.csv")
colnames(pt.df) = paste("PT", colnames(pt.df), sep="-")

# build amalgamated data frame of MT+PT data
amal.df = data.frame()
for(i in 1:length(df$species)) {
  name = df$species[i]
  mt.ref = which(mt.df$Species==name)[1]
  pt.ref = which(pt.df$Species==name)[1]
  if(length(mt.ref)==0 | length(pt.ref)==0) { 
    print(name) 
  } else {
    amal.df = rbind(amal.df, cbind(data.frame(species=name), mt.df[i,2:ncol(mt.df)], pt.df[i,2:ncol(pt.df)]))
  }
}

amal.df
amal.mat = as.matrix(amal.df[,2:ncol(amal.df)])
amal.mat.uniq = unique(amal.mat)

# cluster and visualise
ph = pheatmap(amal.mat.uniq,
              clustering_distance_cols = "binary",
              clustering_distance_rows = "binary", 
              color=c("#000000", "#FFFFFF"), 
              show_rownames = FALSE,
              legend=FALSE,
              treeheight_row=0, cex=0.5)

sf = 2
png("odna-heatmap.png", width=1200*sf, height=1200*sf, res=72*sf)
ph
dev.off()

ph.tree = as.phylo(ph$tree_col)
g.tree.cluster = ggtree(ph.tree, layout='circular') + 
  geom_tiplab(aes(color=substr(label,1,2)), size=2) +
  theme(legend.position="none")


clusters = cutree(ph$tree_col, k=10)
which(clusters==1)

######## oncotree treatment of the loss behaviour
otree = oncotree.fit(1-amal.mat.uniq)
otree.df = data.frame(from=otree$parent$parent, to=otree$parent$child)
otree.g = graph_from_data_frame(otree.df)
V(otree.g)$organelle = "MT"
V(otree.g)$organelle[grepl("P", V(otree.g)$name, fixed=TRUE)] = "PT"

g.tree.oncotree = ggraph(otree.g, layout="dendrogram") + 
  geom_edge_diagonal(edge_width = 1,alpha=0.2) + 
  geom_node_point(aes(color=organelle)) + 
  geom_node_text(aes(label=name, color=organelle),hjust=0,nudge_x=0.8,angle=45,size=2) + theme_void()

sf = 2
png("fig-3.png", width=1200*sf, height=400*sf, res=72*sf)
ggarrange(g.tree.cluster, g.tree.oncotree, nrow=1, widths=c(1,2), labels=c("A", "B"), font.label=list(size=18))
dev.off()

# get the dependent set
start_node = which(V(otree.g)$name == "MT-atp9")
all.paths = all_shortest_paths(otree.g, from = start_node, to = V(otree.g))
# Extract the nodes from the paths
reachable.nodes <- unlist(all.paths$res)
# Remove duplicates and the starting node itself
reachable.nodes <- setdiff(unique(reachable.nodes), start_node)
rn.labels = V(otree.g)$name[reachable_nodes]
pt.set = c()
for(i in 1:length(rn.labels)) {
  tmp = strsplit(rn.labels[i], split="-")[[1]][2]
  pt.set = c(pt.set, tmp)
}

# pull retention indices of these sets
mt.interest = c("ccmb", "rps1", "atp4", "atp1", "rps3", "atp9")
pt.interest = pt.set
mt.df = read.csv("mt-barcode-manual-indices.csv")
pt.df = read.csv("pt-barcode-manual-indices.csv")

mt.df$Subset = 0
mt.df$Subset[mt.df$GeneLabel %in% mt.interest] = 1
mt.df$Subset = factor(mt.df$Subset)

pt.df$Subset = 0
pt.df$Subset[pt.df$GeneLabel %in% pt.interest] = 1
pt.df$Subset = factor(pt.df$Subset)

g.mt = ggplot(mt.df, aes(x=reorder(GeneLabel, -Index), y=Index, fill=Subset)) +
  geom_bar(stat="identity") + theme_light() +
  theme(axis.text.x = element_text(angle=90)) + xlab("") +
  scale_fill_manual(values = c("grey", "black")) 
g.pt = ggplot(pt.df, aes(x=reorder(GeneLabel, -Index), y=Index, fill=Subset)) +
  geom_bar(stat="identity") + theme_light() +
  theme(axis.text.x = element_text(size=3, angle=90)) + xlab("") +
  scale_fill_manual(values = c("grey", "black"))

sf = 2
png("fig-3-supp.png", width=800*sf, height=400*sf, res=72*sf)
ggarrange(g.mt, g.pt, nrow=2, labels=c("MT", "PT"))
dev.off()

mt.n.1 = which(colnames(amal.mat.uniq)=="MT-atp9")
pt.n.1 = which(colnames(amal.mat.uniq)=="PT-ndhf")

amal.mat.uniq[,c(mt.n.1, pt.n.1)]
####### skeletons algorithm
write.table(amal.mat.uniq, "amal-mat-uniq.csv", sep=",", row.names=FALSE, col.names=FALSE)
system("python3 algo3.py amal-mat-uniq.csv")

fname = "amal-mat-uniq.csv-outs2.csv"
df = read.csv(fname, colClasses=c("character"))

g = simplify(graph.data.frame(df))
degs = degree(g)
outdegs = degree(g, mode="out")
n.leaf = length(degs[degs==1])
e.branch = sum(outdegs[outdegs!=0]-1)
tstr = paste(c(e.branch, " excess branches, ", n.leaf, " leaves"), collapse="")
if(length(V(g)) < 50) { vsize = 3 } else { vsize =2 }
#vsize = 2
xoff = 1
g.skeleton = ggraph(g, layout="dendrogram", circular=FALSE) + 
  geom_edge_diagonal() + geom_node_point() + 
  #geom_node_text(aes(label=name),hjust=0,nudge_x=xoff,angle=45,size=vsize) + 
  ggtitle(tstr) + theme_void()
g.skeleton

###############
########## bug hunting
# laminaria digitatus has a weirdly low PT count
mydf3[mydf3$clade=="Phaeophyceae",]

# the dataset we're using doesn't match the intersection from our project
mt.set = mt.df$`MT-Species`
pt.set = pt.df$`PT-Species`

both.set = intersect(mt.set, pt.set)

length(intersect(both.set, mydf3$label))

grepl("laminaria", both.set)

both.tree.labels = c(both.tree$tip.label, both.tree$node.label)
rhodo.ref = which(both.tree.labels == "Rhodophyta")
ggtree(both.tree, layout='circular', branch.length='none', alpha=1) %<+% mydf3 +
  geom_fruit(geom=geom_col, pwidth=1, alpha=0.2,
             mapping = aes(x=pt, y=label), 
             size=0.2, color="#FF0000") +
  geom_fruit(geom=geom_col, pwidth=1, offset=-1, alpha=0.2,
             mapping = aes(x=mt, y=label), 
             size=0.2, color="#0000FF") +
  geom_cladelabel(node=rhodo.ref, label="Rhodophyta")# +
#geom_hilight(mapping=aes(subset = node %in% which(tree.labs %in% clades.interest), fill=nodelab(both.tree, node)),
#              type = "gradient", gradient.direction = 'rt',
#             alpha = .2) #+
