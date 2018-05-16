heatmap_geneset <- function(data, species="Homo sapiens", category="H", subset, ...){
	sets <- msigdbr(species=species, category=category)
	sub_set <- sets[sets$gs_name==subset,]$gene_symbol
	x <- data[row.names(data) %in% sub_set,]
	x <- data.frame(t(apply(x, 1, rescale)))
	heatmap_GSEA(data=x, ...)
	}



	
