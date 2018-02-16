#cat z_serial_comparison.dat | ./d.make_tabular -s '#' | grep -v "^---" > z_serial_benchmarks.tbl

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)


#######################
#SERIAL ALGORITHMS

dfserial <- read.table('z_serial_benchmarks.tbl', header=TRUE, sep="|")

dfserial <- dfserial %>% select(-Rep,-host,-hash)                          %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm)) 

widedf <- dfserial

dfserial <- melt(dfserial, id.vars=c("Algorithm", "Size", "Steps")) %>%
      mutate(variable=gsub('Step(.)','Step \\1:',gsub('_',' ',variable)))

dfserial <- ddply(dfserial, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
dfserial$variable <- factor(dfserial$variable, levels = rev(sort(as.character(unique(dfserial$variable)))))

temp <- dfserial %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total calculation time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + theme_grey(base_size=6) + theme(legend.position="bottom", legend.margin=margin(t = -0.3, unit='cm')) + guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))
ggsave(filename="img_serial_comparison.pdf", plot=p, units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)









###########################
#SIMPLE PARALLEL ALGORITHMS


df_p <- read.table('z_simple_parallel_comet.tbl', header=TRUE, sep="|")

df_p <- df_p %>% select(-Rep,-host,-hash,-Random_seed)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

widedf <- df_p

df_p <- melt(df_p, id.vars=c("Algorithm", "Size", "Steps")) %>%
      mutate(variable=gsub('Step(.)','Step \\1:',gsub('_',' ',variable)))

df_p <- ddply(df_p, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df_p$variable <- factor(df_p$variable, levels = rev(sort(as.character(unique(df_p$variable)))))

#ggplot(df_p, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)")
temp <- df_p %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total calculation time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + theme_grey(base_size=6) + theme(legend.position="bottom", legend.margin=margin(t = -0.3, unit='cm')) + guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))
ggsave(filename="img_simple_parallel_comparison.pdf", plot=p, units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)











################################
#*+PI.exe : Improved Parallelism
################################

df_pi <- read.table('z_parallel_improved_comet.tbl', header=TRUE, sep="|")

df_pi <- df_pi %>% select(-Rep,-host,-hash,-Random_seed)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

widedf <- df_pi

df_pi <- melt(df_pi, id.vars=c("Algorithm", "Size", "Steps")) %>%
      mutate(variable=gsub('Step(.)','Step \\1:',gsub('_',' ',variable)))

df_pi <- ddply(df_pi, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df_pi$variable <- factor(df_pi$variable, levels = rev(sort(as.character(unique(df_pi$variable)))))

#ggplot(df_pi, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)")
temp <- df_pi %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total calculation time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + theme_grey(base_size=6) + theme(legend.position="bottom", legend.margin=margin(t = -0.3, unit='cm')) + guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))
ggsave(filename="img_parallel_improved_comparison.pdf", plot=p, units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)









#################################
#*+PQ.exe : Separated parallelism
#################################

df_pq <- read.table('z_parallel_sep_thread_comet.tbl', header=TRUE, sep="|")

df_pq <- df_pq %>% select(-Rep,-host,-hash,-Random_seed)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

widedf <- df_pq

df_pq <- melt(df_pq, id.vars=c("Algorithm", "Size", "Steps")) %>%
      mutate(variable=gsub('Step(.)','Step \\1:',gsub('_',' ',variable)))

df_pq <- ddply(df_pq, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df_pq$variable <- factor(df_pq$variable, levels = rev(sort(as.character(unique(df_pq$variable)))))

#ggplot(df_pq, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + scale_fill_manual(values=c('#00BFC4'), guide = guide_legend()) + theme(legend.position="bottom")
temp <- df_pq %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total calculation time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + scale_fill_manual(values=c('#00BFC4'), guide = guide_legend()) + theme_grey(base_size=6) + theme(legend.position="bottom", legend.margin=margin(t = -0.3, unit='cm')) + guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))
ggsave(filename="img_parallel_sep_thread_comet.pdf", plot=p, units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)







#################################
#*+GPU.exe : GPU
#################################

df_gpu <- read.table('z_gpu_summitdev.tbl', header=TRUE, sep="|")

df_gpu <- df_gpu %>% select(-Rep,-host,-hash,-Random_seed)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

widedf <- df_gpu

df_gpu <- melt(df_gpu, id.vars=c("Algorithm", "Size", "Steps")) %>%
      mutate(variable=gsub('Step(.)','Step \\1:',gsub('_',' ',variable)))

df_gpu <- ddply(df_gpu, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df_gpu$variable <- factor(df_gpu$variable, levels = rev(sort(as.character(unique(df_gpu$variable)))))

#ggplot(df_gpu, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + scale_fill_manual(values=c('#00BFC4'), guide = guide_legend()) + theme(legend.position="bottom")
temp <- df_gpu %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total calculation time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + scale_fill_manual(values=c('#00BFC4'), guide = guide_legend()) + theme_grey(base_size=6) + theme(legend.position="bottom", legend.margin=margin(t = -0.3, unit='cm')) + guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))
ggsave(filename="img_gpu_summitdev.pdf", plot=p, units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)





#################################
#*+GPU.exe : GPU Scaling
#################################

df_gpu_scaling <- read.table('z_gpu_scaling_summitdev.tbl', header=TRUE, sep="|")

df_gpu_scaling <- df_gpu_scaling %>% select(-Rep,-host,-hash,-Random_seed)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

df_gpu_scaling <- melt(df_gpu_scaling, id.vars=c("Algorithm", "Size", "Steps")) %>%
      mutate(variable=gsub('Step(.)','Step \\1:',gsub('_',' ',variable)))

df_gpu_scaling <- ddply(df_gpu_scaling, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df_gpu_scaling$variable <- factor(df_gpu_scaling$variable, levels = rev(sort(as.character(unique(df_gpu_scaling$variable)))))

#ggplot(df_gpu_scaling, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)") + scale_fill_manual(values=c('#00BFC4'), guide = guide_legend()) + theme(legend.position="bottom")
temp <- df_gpu_scaling %>% filter(variable=="Overall")

p<-ggplot(temp, aes(x=Size**2, y=mean, color=Algorithm)) + 
     #geom_rect(aes(xmin = 0,    xmax = 400, ymin = min(temp$mean), ymax = Inf), fill = "pink", alpha = 0.01, color=NA) +
     #geom_rect(aes(xmin = 400,  xmax = 1000, ymin = min(temp$mean), ymax = Inf), fill = "pink", alpha = 0.01, color=NA) +
     #geom_rect(aes(xmin = 1000, xmax = 2500, ymin = min(temp$mean), ymax = Inf), fill = "pink", alpha = 0.01, color=NA) +
     geom_vline(xintercept=400**2)                                              +
     geom_vline(xintercept=1000**2)                                             +
     geom_vline(xintercept=2500**2)                                             +
     geom_text(aes(x=200**2,  y=80, label="A"), color="black")                  +
     geom_text(aes(x=750**2,  y=80, label="B"), color="black")                  +
     geom_text(aes(x=1750**2, y=80, label="C"), color="black")                  +
     geom_text(aes(x=6500**2, y=80, label="D"), color="black")                  +
     geom_point()                                                               + 
     geom_line()                                                                + 
     scale_x_log10()                                                            + 
     scale_y_log10()                                                            + 
     annotation_logticks()                                                      + 
     theme(axis.text.x=element_text(angle=45, hjust=1))                         + 
     xlab("Cells")                                                              + 
     ylab("Wall-time (s)")                                                      + 
     scale_color_manual(values=c('#00BFC4'), guide = guide_legend())            + 
     theme_grey(base_size=6)                                                    + 
     theme(legend.position="bottom", legend.margin=margin(t = -0.3, unit='cm')) + 
     guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))                   

ggsave(filename="img_gpu_scaling_summitdev.pdf", plot=p, units="in", width=3, height=1.5)

lmdat <- temp %>% filter(0<=Size & Size<=400)
lm(log10(lmdat$mean) ~ log10(lmdat$Size**2))

lmdat <- temp %>% filter(400<=Size & Size<=1000)
lm(log10(lmdat$mean) ~ log10(lmdat$Size**2))

lmdat <- temp %>% filter(1000<=Size & Size<=2500)
lm(log10(lmdat$mean) ~ log10(lmdat$Size**2))

lmdat <- temp %>% filter(2500<=Size)
lm(log10(lmdat$mean) ~ log10(lmdat$Size**2))






############################
#Overall timing comparison



dfc <- rbind(dfserial,df_p,df_pi,df_pq,df_gpu)

temp <- dfc %>% filter(variable=="Overall" & (Size==1000 | Size==10000)) %>% arrange(desc(mean))

temp$Algorithm      <- factor(temp$Algorithm, levels = c(" RB+GPU ", " RB+PQ ", " RB+PI  ", " B&W+PI ", " RB+P  ",  " B&W+P ",  " RB  ", " B&W " ))
temp$Implementation <- temp$Algorithm
temp$Algorithm      <- gsub(' ','',gsub('\\+.*','',temp$Algorithm))



p<-ggplot(temp, aes(x=Implementation, y=mean, fill=Algorithm)) + 
   facet_wrap(~Size) + 
   geom_bar(stat="identity", position="dodge") +
   theme(axis.text.x=element_text(angle=45, hjust=1)) + 
   facet_wrap(~Size, scales="free_x") + 
   coord_flip() + 
   xlab("") + 
   ylab("Wall-time (s)") +
   # scale_fill_manual(values=c('#00BFC4'), guide = guide_legend()) + 
   theme_grey(base_size=6) +
   theme(legend.position="bottom", legend.margin=margin(t = -0.3, unit='cm')) + 
   guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))
ggsave(filename="img_timing_compare.pdf", plot=p, units="in", width=3, height=1.5)
