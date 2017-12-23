#cat z_serial_comparison.dat | ./d.make_tabular -s '#' | grep -v "^---" > z_serial_benchmarks.tbl

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)





df <- read.table('z_serial_benchmarks.tbl', header=TRUE, sep="|")

df <- df %>% select(-Rep,-host,-hash)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

df <- melt(df, id.vars=c("Algorithm", "Size", "Steps"))

df <- ddply(df, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable) %>% mutate(mean=mean/1e6, sd=sd/1e6)
df$variable <- factor(df$variable, levels = sort(as.character(unique(df$variable))))
#df$Size     <- factor(df$Size)

df %>% group_by(variable,Size,Steps) %>% summarise(meandiff=max(mean)-min(mean), sddiff=max(sd)-min(sd), pmean=(max(mean)-min(mean))/max(mean))

ggplot(df, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_y")

temp <- df %>% filter(variable=='Overall')
ggplot(temp, aes(x=Size, y=mean, color=Algorithm)) + geom_point() + geom_smooth(method = "lm", se=FALSE) + theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_y_log10()




df %>% group_by(variable,Size,Steps) %>% summarise(meandiff=max(mean)-min(mean), sddiff=max(sd)-min(sd), pmean=(max(mean)-min(mean))/max(mean))
