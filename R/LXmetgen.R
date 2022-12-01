
LXmetgen <- function(gene_file,meta_file,group1,group2){

  # list all the packages that have been installed
  all_packages <- data.frame(installed.packages())

  # To judge whether a package was installed. If not, it will be installed.
  pack <- data.frame(c("pak","devtools","BiocManager","plyr","eoffice","openxlsx","dplyr","psych","ggplot2",
            "ggrepel","VennDiagram","ggvenn","RColorBrewer","ggthemes","rticles","httr","magrittr",
            "roxygen2","XML","RCurl","curl","stringr","momr","patchwork","ggpubr","scales","conflicted") )

  pack$type <- pack[,1] %in% all_packages$Package

  for (i in 1:nrow(pack)){
    if(pack[i,2]==FALSE)
      install.packages(pack[i,1])
    }
  rm(i)


  # 批量library
  packages <- as.character(pack[,1])

  for(i in packages){
    library(i, character.only = T)
  }
  rm(i)



  #-------- Install the development version of "tidyverse" from GitHub
  if("tidyverse" %in% all_packages$Package==FALSE)
      pak::pak("tidyverse/tidyverse")
     library(tidyverse)


  BiocManager_pack <- data.frame(c("DOSE","clusterProfiler","do","enrichplot",
                        "pathview","BiocParallel","org.Hs.eg.db"))

  BiocManager_pack$type <- BiocManager_pack[,1] %in% all_packages$Package

  for (i in 1:nrow(BiocManager_pack)){
    if(BiocManager_pack[i,2]==FALSE)
      BiocManager::install(BiocManager_pack[i,1])
     }

  # 批量library
  Bio_packages <- as.character(BiocManager_pack[,1])
  for(i in Bio_packages){
    library(i, character.only = T)
  }
  rm(i)


#------------------Gene enriched pathways analysis-----------------------------#

if(dir.exists("temporary files")==FALSE)
    dir.create("temporary files")

if(dir.exists("analysis results")==FALSE)
    dir.create("analysis results")

group <- data.frame(group1,group2)

write.xlsx(group,"temporary files/group.xlsx")

group <- paste("(",group1,"VS",group2,")")

gene_df_0 <- read.xlsx(gene_file)

gene_df <- data.frame(distinct(gene_df_0, gene_df_0[,1], .keep_all = TRUE))

gene_id <- toupper(gene_df[,1])

g_frame <- data.frame(gene_id)

gene_ENTREZID<-bitr(gene_id,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
gene_ENTREZID <- na.omit(gene_ENTREZID)

kegg_gene_df <- enrichKEGG(gene_ENTREZID$ENTREZID, organism = 'hsa',
                           keyType = 'kegg', pvalueCutoff = 0.05,
                           pAdjustMethod = 'BH', minGSSize = 3,
                           maxGSSize = 3500, qvalueCutoff = 0.2,
                           use_internal_data = FALSE)

kegg_gene_result <- kegg_gene_df@result

n_path <- nrow(kegg_gene_result) %>% as.numeric()

path_n <- case_when(n_path>30 ~30,
                    TRUE ~n_path)


title_gene_text <- case_when(n_path>30 ~paste("TOP 30 gene enriched pathways","(",group1,"VS",group2,")"),
                             TRUE ~ paste("Gene enriched pathways","(",group1,"VS",group2,")")
                            )

title_meta_text <- paste("Metabolites KEGG pathways","(",group1,"VS",group2,")")


gene_title_size <- case_when(path_n>=30 ~10,
                             path_n>=20 ~11,
                             TRUE ~12)

gene_xy_size <- case_when(path_n>=30 ~9,
                         path_n>=20 ~10,
                         TRUE ~11)

mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =gene_title_size),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(plot.title = element_text(hjust = 0.5))

xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=10,angle =0,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=gene_xy_size))+
  theme(legend.text=element_text(face="bold",color="black",size=gene_xy_size))

gene_point <- ggplot(kegg_gene_df,showCategory=path_n)+
              geom_point(aes(x=GeneRatio,y=fct_reorder(Description,Count),
              color=-log10(pvalue),size=Count))+
              scale_color_gradient2(midpoint = 1,low = "blue",mid = "#ffcc00",high ="#dc143c")+
              labs(x = 'GeneRatio', y = '',title=title_gene_text)+
             #theme(plot.title = element_text(hjust = 0.5))+
             #theme(axis.title =  element_text(size=10,face = "bold"),
             #axis.text.x = element_text(size=10))+
             mytheme+xytheme+
             theme(panel.grid =element_line(colour="#dcdcdc",size=0.2,linetype = "dashed"))+
             theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

gene_point

ggsave("analysis results/Gene_enriched_pathways.png",gene_point,width=1200, height =1000, dpi=150,units = "px")

df1 <- kegg_gene_df@result

KEGG_T_Symbol <- setReadable(kegg_gene_df, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
df2 <- KEGG_T_Symbol@result

kegg_gene_pathway <- cbind(df1[,1:8],df2[,8],df1[,9])

colnames(kegg_gene_pathway) <- c("ID","Description","GeneRatio","BgRatio",
                                 "pvalue","p.adjust","qvalue","geneID","geneSymbol","Count")

write.xlsx(kegg_gene_pathway,file = "temporary files/gene kegg pathway.xlsx")

gene_path_0 <- read.xlsx("temporary files/gene kegg pathway.xlsx",startRow = 1,rowNames = T)
gene_path_result <- cbind(gene_path_0$Description,gene_path_0$pvalue)
colnames(gene_path_result) <- c("pathways","gene_pvalue")
write.xlsx(data.frame(gene_path_result),"analysis results/gene_enriched_pathway_data.xlsx")


#----------------Metabolite enriched pathways analysis-------------------------#

# install the dependency packages using a function [metanr_packages()]
metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz",
                 "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph",
                 "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR",
                 "fgsea", "crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs,force=TRUE,quietly = TRUE)
    print(c(new_pkgs, " packages added..."))
   }

  if((length(new_pkgs)<1)){
    print("No new packages added...")
   }
 }


metanr_packages()

# To judge whether a package was installed. If not, it will be installed.
if('MetaboAnalystR' %in% all_packages$Package==FALSE)
  install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)

  library(MetaboAnalystR)


# reading the data file
meta_set <- read.xlsx(meta_file)

# To judge the type of the data, "HMDB" or "name"
if(substr(meta_set[1,1],1,4)=="HMDB")
  inputType = "hmdb" else
    inputType = "name"

# translate the first column as a character list
m_set <- t(meta_set[,1]) %>% as.character()

# perform the packages below to get the mapping metabolites from the KEGG database.
m_obj<-InitDataObjects("conc", "msetora", FALSE)
m_obj<-Setup.MapData(m_obj, m_set)


url <- "https://www.metaboanalyst.ca/resources/libs/compound_db.qs"
if(file.exists("compound_db.qs")==FALSE)
curl_download(url, "compound_db.qs")

m_obj<-CrossReferencing(m_obj, inputType)
m_obj<-CreateMappingResultTable(m_obj)
m_Filter<-SetMetabolomeFilter(m_obj, F)

# screen out the mapping result from the m_Filter list.
mapping_result_all <- m_Filter$dataSet$map.table %>% data.frame()

# screen out the mapping result without the KEGG ID "NA".
mapping_result_all$KEGG[which(mapping_result_all$KEGG=="NA")]=NA
mapping_result <- dplyr::filter(mapping_result_all,!is.na(KEGG))

write.xlsx(mapping_result_all,"analysis results/metabolite_mapping_data_all.xlsx")
write.xlsx(mapping_result,"analysis results/metabolite_mapping_data.xlsx")

m_lib<-SetKEGG.PathLib(m_obj, "hsa", "current")
m_filt<-SetMetabolomeFilter(m_lib, F);
m_scor<-CalculateOraScore(m_filt, "rbc", "hyperg")

path_result <- read.csv("pathway_results.csv",header=T)

file.remove("pathway_results.csv")

colnames(path_result) <- c("path_id","Total","Expected","Hits","p value","-log10(p)","Holm p","FDR","Impact")


if(file.exists("analysis results/kegg_hsa_pathways.xlsx")==FALSE){
  website <- "https://rest.kegg.jp/list/pathway/hsa"

  t <- url(website,encoding ="UTF-8")  # library(XML)
  d<-scan(t, what=character()) %>% paste0(collapse = " ")
  d_c <- gsub("- Homo sapiens (human)",";",d,fixed = TRUE)
  dw <- strsplit(d_c, "[;]") %>% data.frame()
  dw[,1] <- trimws(dw[,1])
  dw$path_id <- substring(dw[,1],6,14)
  dw$path_name <- substring(dw[,1], 15, str_length(dw[,1]))  %>% trimws()
  hsa_paths <- dw[-1]
  write.xlsx(hsa_paths,"analysis results/kegg_hsa_pathways.xlsx")
 }

path_ls <- read.xlsx("analysis results/kegg_hsa_pathways.xlsx",colNames=T) %>% data.frame()

path_ls$path_id <- trimws(as.character(path_ls$path_id))

path_result$path_id <- trimws(as.character(path_result$path_id))

path_df <- left_join(path_result,path_ls,by="path_id")

meta_path_data <-cbind(path_df[1],path_df[10],path_df[,2:9])

meta_path_data[,7] <- -log10(meta_path_data[,6])

colnames(meta_path_data) <- c("path_id","path_name","Total","Expected","Hits","Pvalue","-log10(p)","Holm p","FDR","Impact")

write.xlsx(meta_path_data,"analysis results/metabolite_enriched_pathway_data.xlsx",rownames=T,colnames=T)


title_size <- case_when(nrow(meta_path_data)>30 ~13,
                        nrow(meta_path_data)>20 ~12,
                        TRUE ~11)

xy_size <- case_when(nrow(meta_path_data)>30 ~8,
                     nrow(meta_path_data)>20 ~9,
                     TRUE ~10)

my_theme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"))+
  theme(plot.title = element_text(hjust = 0.5))


xy_theme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size,angle =0,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=xy_size))+
  theme(legend.text=element_text(face="bold",color="black",size=10))


x_ran <- max(meta_path_data$Impact)*1.05
y_ran <- nrow(meta_path_data)

meta_path <- ggplot(meta_path_data,
                    aes(x=Impact,
                        y=fct_reorder(path_name,-log10(Pvalue))))+
  geom_point(aes(color=-log10(Pvalue),size=Impact))+
  scale_color_gradient2(midpoint = 0,low = "#33ffff",mid = "#ffff00",high = "#ff0000")+
  labs(x = 'Pathway Impact', y = '',title="Overview of Metabolism Pathway Anlysis")+
  my_theme+xy_theme+
  theme(panel.grid =element_line(colour="#dcdcdc",size=0.2,linetype = "dashed"))+
  coord_cartesian(xlim=c(0, x_ran),expand = T)+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

meta_path

ggsave("analysis results/Metabolite_enriched_Pathways.png", meta_path, width=1200, height =1000, dpi=150,units = "px")

all_files <- dir()
qs <- grep("*.qs",all_files)
rds <- grep("*.rds",all_files)
file.remove(all_files[rds],"name_map.csv")


#---------Jiont analysis of the gene and metabolite enriched pathways---------#

meta_path_0 <- meta_path_data[,-1]

colnames(meta_path_0) <- c("pathways","Total","Expected","Hits",
                           "pvalue","-log10(p)","Holm.adjust","FDR","Impact")

meta_path_result <- cbind(meta_path_0$pathways,meta_path_0$pvalue)
colnames(meta_path_result) <- c("pathways","meta_pvalue")

set_gene <- data.frame(gene_path_result[,1])
set_meta <- data.frame(meta_path_result[,1])

venn <- venn.diagram(
  c(set_gene,set_meta),show_percentage=TRUE,filename=NULL,
  main="Venn graphics",main.pos = c(0.5, 1.05), main.fontface = "bold",
  main.col = "black", main.cex = 2, main.just = c(0.5, 1),
  sub.pos = c(0.5,1.05), sub.fontface = "plain", sub.fontfamily = "serif",
  sub.col = "black", sub.cex = 1, sub.just =c(0.5, 1),
  category.names = c(group1,group2),
  col=c("turquoise", 'sandybrown'),
  lwd=1,
  #lty="dashed",
  fill = c("turquoise",'sandybrown'),
  cat.col = c("#440154ff", '#21908dff'),
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 2,
  cat.default.pos = "text",
  cat.pos = c(-20, 20),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans"
)

#grid.draw(venn)

file.remove(list.files(pattern = "*log"))

gene_meta_Venn <- inner_join(data.frame(gene_path_result),data.frame(meta_path_result),
                             by="pathways")

#write.xlsx(gene_meta_Venn,"new_data_gene_meta_Venn.xlsx",rowNames=TRUE)
colnames(gene_meta_Venn) <- c("Pathway","P_gene","P_meta")

gene_path <- gene_meta_Venn[c("Pathway","P_gene")]
gene_P <- as.numeric(gene_path$P_gene)
class(gene_P)
gene_path$P_gene <- -log2(gene_P)

gene_path$types <- rep("genes",nrow(gene_path))
colnames(gene_path) <- c("Pathway","-log2(Pvalue)","types")

meta_path <- gene_meta_Venn[c("Pathway","P_meta")]
meta_P <- as.numeric(meta_path$P_meta)
meta_path$P_meta <- -log2(meta_P)
meta_path$types <- rep("metabolites",nrow(meta_path))
colnames(meta_path) <- c("Pathway","-log2(Pvalue)","types")

gene_meta_path <- bind_rows(gene_path,meta_path)
colnames(gene_meta_path) <- c("Pathways","minus_log2_Pvalue","types")

write.xlsx(gene_meta_path,"temporary files/gene_meta_path.xlsx")

height_y <- max(gene_meta_path$minus_log2_Pvalue)+1

f0 <- ggplot(gene_meta_path, aes(x = Pathways, y = minus_log2_Pvalue,fill=types))
f1 <- f0+
  scale_y_continuous(expand = c(0, 0),limits = c(0, height_y))+
  geom_bar(position = "dodge",stat = "identity")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(colour = "black", face="bold",size=12)) +
  labs(x="",y = ("-log2(Pvalue)"),title = paste('Pathways analysis',group))

f1

log05 <- -log2(0.05)
log01 <- -log2(0.01)

line1 <- geom_hline(yintercept = c(log05),
                    linewidth = 0.6,
                    color = "blue",
                    lty = "dashed")
line2 <- geom_hline(yintercept = c(log01),
                    linewidth = 0.6,
                    color = "red",
                    lty = "dashed")

y1 <- geom_text(x=nrow(gene_meta_path)/2-2.5,y=log05+0.8,label = c("p<0.05"),
                size=4,color="blue",fontface="italic")
y2 <- geom_text(x=nrow(gene_meta_path)/2-2.5,y=log01+0.8,label = c("p<0.01"),
                size=4,color="blue",fontface="italic")

f2 <- f1+line1+line2+y1+y2
f2


nrow_path <- nrow(gene_meta_path)/2

joint_title_size <- case_when(nrow_path>=30 ~12,
                              nrow_path>=20 ~12,
                              TRUE ~12)

joint_x_size <- case_when(nrow_path>=20 ~8,
                              nrow_path>=10 ~9,
                              TRUE ~10)

joint_y_size <- case_when(nrow_path>=20 ~11,
                          nrow_path>=10 ~11,
                          TRUE ~11)

joint_legend_size <- case_when(nrow_path>=30 ~11,
                          nrow_path>=20 ~11,
                          TRUE ~11)
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="black",face="bold",size =joint_title_size),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=2), "cm")
       )+
  theme(plot.title = element_text(hjust = 0.5))


xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=joint_x_size,angle =45,hjust=1))+
  theme(axis.text.y = element_text(face="bold",color="black",size=joint_y_size))

legend_theme <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = joint_legend_size, face = "bold"),
  legend.direction = "vertical",
  #legend.position = c(0.5,0.9),
  legend.background = element_blank()
)

f3 <- f2+mytheme+xytheme

f3

ggsave("analysis results/The gene_metabolite_Joint_pathways 01.png",f3,width=1200, height =1000, dpi=150,units = "px")


f4 <- f3+theme(legend.direction = "horizontal",legend.position = c(0.5,0.9))+labs(fill="")

ggsave("analysis results/The gene_metabolite_Joint_pathways 02.png",f4,width=1200, height =1000, dpi=150,units = "px")

print("--------------------------------------------------------------")
print("The results can be found in the folder of <analysis results>")

f4

}



