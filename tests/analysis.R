#cat z_serial_comparison.dat | ./d.make_tabular -s '#' | grep -v "^---" > z_serial_benchmarks.tbl

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)

df <- read.table('z_serial_benchmarks.tbl', header=TRUE, sep="|")

df <- df %>% select(-Rep,-host,-hash)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

widedf <- df

df <- melt(df, id.vars=c("Algorithm", "Size", "Steps"))

df <- ddply(df, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df$variable <- factor(df$variable, levels = sort(as.character(unique(df$variable))))
#df$Size     <- factor(df$Size)






temp <- df %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total_calculation_time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)")
ggsave(filename="img_serial_comparison.pdf", plot=p+guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))+theme_grey(base_size=6), units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)










df <- read.table('z_simple_parallel_comet.tbl', header=TRUE, sep="|")

df <- df %>% select(-Rep,-host,-hash,-Random_seed)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

widedf <- df

df <- melt(df, id.vars=c("Algorithm", "Size", "Steps"))

df <- ddply(df, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df$variable <- factor(df$variable, levels = sort(as.character(unique(df$variable))))
#df$Size     <- factor(df$Size)


ggplot(df, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)")
temp <- df %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total_calculation_time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)")
ggsave(filename="img_simple_parallel_comparison.pdf", plot=p+guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))+theme_grey(base_size=6), units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)











################################
#*+PI.exe : Improved Parallelism
################################

df <- read.table('z_parallel_improved_comet.tbl', header=TRUE, sep="|")

df <- df %>% select(-Rep,-host,-hash,-Random_seed)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

widedf <- df

df <- melt(df, id.vars=c("Algorithm", "Size", "Steps"))

df <- ddply(df, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df$variable <- factor(df$variable, levels = sort(as.character(unique(df$variable))))
#df$Size     <- factor(df$Size)


ggplot(df, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)")
temp <- df %>% filter(Size==1000 | Size==10000) %>% filter(variable!="Total_calculation_time")
p<-ggplot(temp, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_x") + coord_flip() + xlab("") + ylab("Wall-time (s)")
ggsave(filename="img_parallel_improved_comparison.pdf", plot=p+guides(fill = guide_legend(keywidth=0.4, keyheight=0.4))+theme_grey(base_size=6), units="in", width=3, height=1.5)

#Percentage of time taken by erosion
erosion_per<-widedf %>% select(Algorithm,Steps,Size,Overall,Step7_Erosion) %>% mutate(per=Step7_Erosion/Overall)
mean(erosion_per$per)













































df %>% group_by(variable,Size,Steps) %>% summarise(meandiff=max(mean)-min(mean), sddiff=max(sd)-min(sd), pmean=(max(mean)-min(mean))/max(mean))

ggplot(df, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_y")





temp <- df %>% filter(variable=='Overall')
ggplot(temp, aes(x=Size, y=mean, color=Algorithm)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_log10() + scale_x_log10()

temp <- df
temp$Size <- factor(temp$Size)
ggplot(temp, aes(x=variable, y=mean, fill=Size)) + geom_bar(stat="identity", position="dodge") + scale_y_log10() + facet_wrap(~Algorithm)


df %>% group_by(variable,Size,Steps) %>% summarise(meandiff=max(mean)-min(mean), sddiff=max(sd)-min(sd), pmean=(max(mean)-min(mean))/max(mean))
