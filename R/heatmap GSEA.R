heatmap_GSEA <- function(data, dendro="none", limits=c(min(heatmap_data$expr), max(heatmap_data$expr)), dendro_thick=0.25, dendro_prop_x=0.2, dendro_prop_y=0.2,x_axis_text_angle=45, x_axis_font_size=6,  y_axis_font_size=5, colorbar_height=10, title=NULL, title_size=16, x_y_ratio="square_tile", show_legend=TRUE, dist_method="euclidean", hclust_method="complete") {
	x_y_ratio <- ifelse(x_y_ratio=="square_tile", dim(data)[1]/dim(data)[2], x_y_ratio)
	if (dendro=="none"){
		data$gene <- rownames(data)
	
		heatmap_data <- data %>% 
    	reshape2::melt(value.name = "expr", id.vars = c("gene"))
    
    	p <- (ggplot(heatmap_data, 
    			aes(x=variable, y=gene, fill=expr))
    		+geom_tile(col="black", show.legend=show_legend)
    		+scale_fill_gsea(limits=limits,oob=squish)	
    		+scale_x_discrete(
    			position="bottom",
    			expand = c(0, -0.075))
    		+scale_y_discrete(
    			expand = c(0, -0.075), 
    			position="right") + labs(x = "", y = "") 
    		+theme_bw()
    		+theme(
    			axis.text.x = element_text(size = x_axis_font_size, hjust = 1, angle = x_axis_text_angle, color="black"), 
    			axis.text.y = element_text(size=y_axis_font_size, color="black"),
          		# margin: top, right, bottom, and left
          		plot.margin = unit(c(2, 1, 0, 2), "cm"), 
          		panel.grid.minor = element_blank(), 
          		legend.title=element_blank(), 
          		plot.title=element_text(size=title_size, hjust=0.5),
    			    aspect.ratio=x_y_ratio)
          		+guides(fill=guide_colorbar(barwidth=0.5, barheight=colorbar_height))
          		+ggtitle(title))
          	p
    
    } else if(dendro=="both"){
    	
    	x <- as.matrix(data)
		dd.col <- as.dendrogram(hclust(dist(t(x),method=dist_method), method=hclust_method))
		dd.row <- as.dendrogram(hclust(dist(x,method=dist_method), method=hclust_method))
		dx <- dendro_data(dd.row)
		dy <- dendro_data(dd.col)

		segment_data_x <- with(
    		segment(dx), 
    		data.frame(x = y, y = x, xend = yend, yend = xend))
    
		segment_data_y <- with(
    		segment(dy), 
    		data.frame(x = y, y = x, xend = yend, yend = xend))

		gene_pos_table <- with(
    		dx$labels, 
    		data.frame(y_center = x, gene = as.character(label), height = 1))

		sample_pos_table <- with(
    		dy$labels, 
    		data.frame(x_center = x, sample = as.character(label), width = 1))
    
		heatmap_data <- x %>% 
    		reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
    		left_join(gene_pos_table) %>%
    		left_join(sample_pos_table)

		gene_axis_limits <- with(
    		gene_pos_table, 
    		c(min(y_center - 0.5 * height), max(y_center + 0.5 * height)))
    
		sample_axis_limits <- with(
    		sample_pos_table, 
    		c(min(x_center - 0.5 * width), max(x_center + 0.5 * width)))
    		
    	expand_y <- (max(segment_data_x$x)*0.05)*(0.2/dendro_prop_x)
    	expand_x <- (max(segment_data_y$x)*0.05)*(0.2/dendro_prop_y)
    		
 		p <- (ggplot(heatmap_data, 
				aes(x=x_center, y=y_center, fill=expr, height=height, width=width))
			+geom_tile(col="black", show.legend=show_legend)
			+scale_fill_gsea(limits=limits,oob=squish)
			+scale_x_continuous(
				breaks = sample_pos_table$x_center, 
             	labels = sample_pos_table$sample, 
             	expand = c(0, 0), 
              position="bottom")
          	+scale_y_continuous(
             	breaks = gene_pos_table[, "y_center"], 
             	labels = gene_pos_table$gene,
              limits =gene_axis_limits, 
            	expand = c(0, 0), 			
              position="right") 
        	+labs(x = "", y = "") 
        	+theme_bw()
          	+theme(
             	axis.text.x = element_text(size = x_axis_font_size, hjust = 1, angle = x_axis_text_angle, color="black"), 
            	axis.text.y = element_text(size=y_axis_font_size, color="black"),
          		# margin: top, right, bottom, and left
          		plot.margin = unit(c(1, 1, 0, 1), "cm"), 
          		panel.grid.minor = element_blank(), 
          		legend.title=element_blank(), 
          		plot.title=element_text(size=title_size, hjust=0.5), 
          		aspect.ratio=x_y_ratio)
          	+guides(fill=guide_colorbar(barwidth=0.5, barheight=colorbar_height))
          	+ggtitle(title))
    	
    	plt_dendr_x <- (axis_canvas(p, axis="y") 
    		+geom_segment(data=segment_data_x, 
    			aes(x = x, y = y, xend = xend, yend = yend), 
    			size=dendro_thick, 
    			lineend="square") 
    			+scale_x_reverse(expand = c(0, expand_y)) 
    			+scale_y_continuous(
    				breaks = gene_pos_table$y_center, 
              	labels = gene_pos_table$gene, 
                 	limits = gene_axis_limits, 
                 	expand = c(0, 0)) + 
    			labs(x = "", y = "", colour = "", size = "") +
    			theme_void() + 
    			theme(panel.grid.minor = element_blank()))
    
		plt_dendr_y <- (axis_canvas(p, axis="x")
    			+geom_segment(data=segment_data_y, 
    			aes(x = y, y = x, xend = yend, yend = xend),
    			size=dendro_thick,
    			lineend="square") 
    			+scale_x_continuous(
    				breaks = sample_pos_table$x_center, 
                	labels = sample_pos_table$sample, 
                	limits = sample_axis_limits, 
                	expand = c(0, 0), position="top")
                	+labs(x = "", y = "", colour = "", size = "") 
                	+scale_y_continuous(expand = c(0, expand_x))
                	+theme_void()
                	+theme(
                		panel.grid.minor = element_blank(),
                		plot.margin = unit(c(0, 0, 0, 0), "cm")))
    		
    		
    			p %>% insert_xaxis_grob(plt_dendr_y, grid::unit(dendro_prop_x, "null"), position = "top") %>%
  				insert_yaxis_grob(plt_dendr_x, grid::unit(dendro_prop_y, "null"), position = "left") %>%
  				ggdraw()

    	} else if(dendro=="y"){
    	
    		x <- as.matrix(data)
			dd.row <- as.dendrogram(hclust(dist(x, method=dist_method), method=hclust_method))
			dx <- dendro_data(dd.row)

			segment_data_x <- with(
    			segment(dx), 
    			data.frame(x = y, y = x, xend = yend, yend = xend))
    

			gene_pos_table <- with(
    			dx$labels, 
    			data.frame(y_center = x, gene = as.character(label), height = 1))
    
			heatmap_data <- x %>% 
    			reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
    			left_join(gene_pos_table) 

			gene_axis_limits <- with(
    			gene_pos_table, 
    			c(min(y_center - 0.5 * height), max(y_center + 0.5 * height)))
    			
    		expand_y <- max(segment_data_x$x)*0.05
    		
 			p <- (ggplot(heatmap_data, 
					aes(x=sample, y=y_center, fill=expr, height=height))
				+geom_tile(col="black", show.legend=show_legend)
				+scale_fill_gsea(limits=limits,oob=squish)
				+scale_x_discrete(
             		expand = c(0, 0), 
              		position="bottom")
          		+scale_y_continuous(
             		breaks = gene_pos_table[, "y_center"], 
             		labels = gene_pos_table$gene,
              		limits =gene_axis_limits, 
            		expand = c(0, 0), 			
              		position="right") 
        			+labs(x = "", y = "") 
        			+theme_bw()
          			+theme(
             			axis.text.x = element_text(size = x_axis_font_size, hjust = 1, angle = x_axis_text_angle, color="black"), 
            			axis.text.y = element_text(size=y_axis_font_size, color="black"),
          				# margin: top, right, bottom, and left
          				plot.margin = unit(c(2, 1, 0, 1), "cm"), 
          				panel.grid.minor = element_blank(), 
          				legend.title=element_blank(), 
          				plot.title=element_text(size=title_size, hjust=0.5),
          				aspect.ratio=x_y_ratio)
          			+guides(fill=guide_colorbar(barwidth=0.5, barheight=colorbar_height))
          			+ggtitle(title))
    	
    		plt_dendr_x <- (axis_canvas(p, axis="y") 
    			+geom_segment(data=segment_data_x, 
    				aes(x = x, y = y, xend = xend, yend = yend), 
    				size=dendro_thick, 
    				lineend="square") 
    			+scale_x_reverse(expand = c(0, expand_y)) 
    			+scale_y_continuous(
    				breaks = gene_pos_table$y_center, 
              		labels = gene_pos_table$gene, 
                 	limits = gene_axis_limits, 
                 	expand = c(0, 0)) + 
    			labs(x = "", y = "", colour = "", size = "") 
    			+theme_void()
    			+theme(panel.grid.minor = element_blank()))
    
			p %>% insert_yaxis_grob(plt_dendr_x, grid::unit(dendro_prop_y, "null"), position = "left") %>%
  			ggdraw()
    	
    	} else if(dendro=="x"){
    	
    		x <- as.matrix(data)
			dd.col <- as.dendrogram(hclust(dist(t(x), method=dist_method), method=hclust_method))
			dy <- dendro_data(dd.col)
    
			segment_data_y <- with(
    			segment(dy), 
    			data.frame(x = y, y = x, xend = yend, yend = xend))

			sample_pos_table <- with(
    			dy$labels, 
    			data.frame(x_center = x, sample = as.character(label), width = 1))
    
			heatmap_data <- x %>% 
    			reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
    			left_join(sample_pos_table)
    
			sample_axis_limits <- with(
    			sample_pos_table, 
    			c(min(x_center - 0.5 * width), max(x_center + 0.5 * width)))
    		

	    	expand_x <- max(segment_data_y$x)*0.05
    		
    	
    		p <- (ggplot(heatmap_data, 
					aes(x=x_center, y=gene, fill=expr, width=width))
				+geom_tile(col="black", show.legend=show_legend)
				+scale_fill_gsea(limits=limits,oob=squish)
				+scale_x_continuous(
					breaks = sample_pos_table$x_center, 
                    labels = sample_pos_table$sample,
             		expand = c(0, 0), 
              		position="bottom")
          		+scale_y_discrete( 
            		expand = c(0, -0.075), 			
              		position="right") 
        			+labs(x = "", y = "") 
        			+theme_bw()
          			+theme(
             			axis.text.x = element_text(size = x_axis_font_size, hjust = 1, angle = x_axis_text_angle, color="black"), 
            			axis.text.y = element_text(size=y_axis_font_size, color="black"),
          				# margin: top, right, bottom, and left
          				plot.margin = unit(c(1, 1, 0, 2), "cm"), 
          				panel.grid.minor = element_blank(), 
          				legend.title=element_blank(), 
          				plot.title=element_text(size=title_size, hjust=0.5), 
          				aspect.ratio=x_y_ratio)
          			+guides(fill=guide_colorbar(barwidth=0.5, barheight=colorbar_height))
          			+ggtitle(title))
    	
    		plt_dendr_y <- (axis_canvas(p, axis="x")
    			+geom_segment(data=segment_data_y, 
    			aes(x = y, y = x, xend = yend, yend = xend),
    			size=dendro_thick,
    			lineend="square") 
    			+scale_x_continuous(
    				breaks = sample_pos_table$x_center, 
                	labels = sample_pos_table$sample, 
                	limits = sample_axis_limits, 
                	expand = c(0, 0), position="top")
                	+labs(x = "", y = "", colour = "", size = "") 
                	+scale_y_continuous(expand = c(0, expand_x))
                	+theme_void()
                	+theme(
                		panel.grid.minor = element_blank(),
                		plot.margin = unit(c(0, 0, 0, 0), "cm")))
    	
    		p %>% insert_xaxis_grob(plt_dendr_y, grid::unit(dendro_prop_x, "null"), position = "top") %>%
  			ggdraw()
    	}
    }