rm(list=ls())
require(data.table)
require(ggplot2)
require(cowplot)
require(reshape2)
require(RColorBrewer)
require(plyr)
require(doMC)
require(cowplot)
source("~/Dropbox/singlecell-pipeline/fetal_venule/scripts/rnaseq_funcs_V2.R")
setwd('~/Dropbox/singlecell-pipeline/fetal_venule/e12_finalSubmission/iRPCA_NatureSubmission/')

# load data -----

load('../../exp1_nextseq/data/fetal_venule_exp1_Nextseq.dt_goodCells.RData')
load('../../exp1_nextseq/data/fetal_venule_exp1_Nextseq.cell_info.RData')

load('all_r7/all_r7_artl.Rdata')
load('all_r9/all_r9_stalk.RData')
load('all_r9/all_r9_tip.RData')
load('SVv_SVc/all_SVv_SVc_SVc.RData')
load('endocard_refine_r5/cells_endocardRefine.r5.v2.endocard.Rdata')
load('all_r14/all_r14_SVv.RData')
load('all_r11/all_r11_Mes1.RData')
load('all_r11/all_r11_Mes2.RData')
load('all_r14/all_r14_VV1.RData')
load('all_r14/all_r14_VV2.RData')
load('all_r12/all_r12_col1a1.RData')
load('all_r13/all_r13_foxf1a.RData')
load('all_r1/cells_hbb_hi.RData')

cells.use = c(cells.allr7.artl,
              cells.allr9.stalk,
              cells.allr9.tip,
              cells.allSVv_SVc.SVc,
              cells.allr14.SVv,
              cells.endocardRefine.r5.v2.endocard,
              cells.allr14.VV1,
              cells.allr14.VV2,
              cells.allr11.Mes1,
              cells.allr11.Mes2, 
              cells.allr12.col1a1, 
              cells.allr13.foxf1a,
              cells.hbb.hi)

cell.info = fv1_NS.cell.info[cell.name %in% cells.use]
cell.info[cell.name %in% cells.allr7.artl, subtype := 'Arterial']
cell.info[cell.name %in% cells.allr9.stalk, subtype := 'CV1']
cell.info[cell.name %in% cells.allr9.tip, subtype := 'CV2']
cell.info[cell.name %in% cells.allSVv_SVc.SVc, subtype := 'SVc']
cell.info[cell.name %in% cells.allr14.SVv, subtype := 'SVv']
cell.info[cell.name %in% cells.endocardRefine.r5.v2.endocard, subtype := 'Endocardial']
cell.info[cell.name %in% cells.allr14.VV1, subtype := 'Venous.Valve.1']
cell.info[cell.name %in% cells.allr14.VV2, subtype := 'Venous.Valve.2']
cell.info[cell.name %in% cells.allr11.Mes1, subtype := 'Mesenchymal.1']
cell.info[cell.name %in% cells.allr11.Mes2, subtype := 'Mesenchymal.2']
cell.info[cell.name %in% cells.allr12.col1a1, subtype := 'col1a1.excluded']
cell.info[cell.name %in% cells.allr13.foxf1a, subtype := 'foxf1a.excluded']
cell.info[cell.name %in% cells.hbb.hi, subtype := 'rbc.excluded']

cell.info <- cell.info[!is.na(subtype)]
unique(cell.info[,subtype])

tmp1.cell.info <- cell.info[!subtype %in% c('col1a1.excluded','foxf1a.excluded','rbc.excluded','SVc','Arterial','CV1','CV2')]
tmp1.cell.info[, subtype_color := brewer.pal(8, 'Dark2')[factor(subtype)]]

tmp2.cell.info <- cell.info[subtype %like% 'excluded']
tmp2.cell.info[, subtype_color := gg_color_hue(3)[as.factor(subtype)]]

tmp3.cell.info <- cell.info[subtype %in% c('SVc','Arterial','CV1','CV2')]
tmp3.cell.info[subtype=="SVc", subtype_color:=rgb(119, 207, 244, maxColorValue = 255)]
tmp3.cell.info[subtype=="Arterial", subtype_color:=rgb(242, 122, 129, maxColorValue = 255)]
tmp3.cell.info[subtype=="CV1", subtype_color:=rgb(141,113, 252, maxColorValue = 255)]
tmp3.cell.info[subtype=="CV2", subtype_color:=rgb(246,155,235,maxColorValue = 255)]

cell.info <- rbindlist(list(tmp1.cell.info, tmp2.cell.info, tmp3.cell.info), use.names = T, fill = F)
unique(cell.info[,.(subtype, subtype_color)])

subtype2col = as.character(cell.info[order(subtype), subtype_color])
names(subtype2col) = as.character(cell.info[, subtype])

fv1.subtype.info = copy(cell.info)


save(fv1.subtype.info, file = '../data/e12_subtypeInfo.Rdata')
write.csv(fv1.subtype.info, file = '../data/e12_subtypeInfo.csv')

plt=ggplot(fv1.subtype.info[!subtype %like% 'excluded'], aes(subtype, num.genes, fill=subtype_color)) +
  geom_violin() +
  scale_fill_identity(guide = "none") +
  theme(axis.text = element_text(angle=45, hjust=1))
save_plot('../qc/numgenes_subtype_e12.pdf', plt, base_height = 3, base_aspect_ratio = 3)


# Add subtypes to Seurat data
load('../../exp1_nextseq/data/fetal_venule_e12_seurat.Rdata')
seur@meta.data$subtype <- NULL
seur@meta.data$subtype_color <- NULL
metadata <- as.data.frame(fv1.subtype.info)
rownames(metadata) <- metadata$cell.name
seur <- AddMetaData(seur, metadata = metadata)
save(seur, file="../data/e12_seurat.Rdata")
