## Written 2019 by Asa Bjorklund asa.bjorklund@scilifelab.se
# OBS! Cannot install seurat3 on shiny server, need to do plotting with ggplot.


# function for printing grid
print.myplot <- function(x, ...) {
  grid::grid.draw(x)
}



# some extra settings for ggplot
noax <- theme_void() + theme(plot.margin = unit(c(0,1,0,1), "cm"))
noleg <- theme(legend.position="none")
legsize <- theme(legend.title=element_text(size=10), legend.text=element_text(size=8))

scale <- scale_color_gradientn(colors=c("yellow","red","black"))

# define custom colors for some metadata
col.scale <- list()
col.cluster <-  c("blue1", "gray", "red1","green3")
names(col.cluster) <- c("Astrocyte","Cycling","Neuron","VLMC")
col.scale$CelltypeCluster <- scale_colour_manual(values = col.cluster)
col.scale$Type <- scale_colour_manual(values = c("cyan", "magenta"))



plot.tsne <- function(input){
  plotname <- input$color1
  temp <- meta
  cont <- FALSE
  if (plotname %in% c("gene","nFeature_RNA")) {
     cont <- TRUE
  }
  

  if (input$color1 == "gene"){
     gene <- input$gene1
     temp[gene] <- expr[gene,]
     plotname <- gene
  }

  if (input$reduction == "UMAP"){
     x <- "umap_x"
     y <- "umap_y"
  }else{
     x <- "tsne_x"
     y <- "tsne_y"
  }


  if (cont){
     p1 <- ggplot() + geom_jitter(data=temp,aes_(x=as.name(x),y=as.name(y),color=as.name(plotname)),na.rm=TRUE) + noax + scale + legsize
     p2 <- ggplot() + geom_violin(data=temp,aes_(x=as.name("CelltypeCluster"),y=as.name(plotname),fill=as.name("CelltypeCluster") ), scale="width")+ theme_classic() + noleg + scale_fill_manual(values=col.cluster)
     plot <- arrangeGrob(p1,p2,ncol=2,widths=c(4,3))
  }else{
     plot <- ggplot() + geom_jitter(data=temp,aes_(x=as.name(x),y=as.name(y),color=as.name(plotname)),na.rm=TRUE) + noax + legsize
     if (plotname %in% names(col.scale)){
     	plot <- plot + col.scale[[plotname]]
     }

  }

  class(plot) <- c("myplot", class(plot))
  return(plot)
}
  

######################################
# the actual server

shinyServer(  function(input, output,session) {

  plot <- reactive( { plot.tsne(input) })
  output$plot_tsne <- renderPlot( { print(plot()) } )

  output$download <- downloadHandler(
    filename <- "tsnes_grafting.pdf",
    content = function(file) {
      pdf(file,width=10,height=4)              
      print(plot())
      dev.off()
  })  
})



