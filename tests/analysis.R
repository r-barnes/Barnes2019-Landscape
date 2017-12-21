#./d.make_tabular -s '#' /z/out > ./benchmarks.tbl

library(dplyr)
library(reshape2)
library(plyr)
library(ggplot2)

df <- read.table('benchmarks.tbl', header=TRUE, sep="|")

df <- df %>% select(-Prog,-Rep)                                     %>% 
             mutate(Algorithm=gsub('FastScape_','',Algorithm))

#df <- df %>% group_by(Algorithm,Size) %>% summarise_all(.funs = c(mean="mean", sd="sd"))

df <- melt(df, id.vars=c("Algorithm", "Size", "Steps"))

df <- ddply(df, c("Algorithm", "Size", "Steps", "variable"), summarise, mean = mean(value), sd = sd(value)) %>% arrange(Algorithm,variable)
df$variable <- factor(df$variable, levels = sort(as.character(unique(df$variable))))

ggplot(df, aes(x=variable, y=mean, fill=Algorithm)) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + facet_wrap(~Size, scales="free_y")