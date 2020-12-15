## Written 2018 by Asa Bjorklund asa.bjorklund@scilifelab.se

# function for printing grid
print.myplot <- function(x, ...) {
  grid::grid.draw(x)
}



# some extra settings for ggplot

noax <- theme_void() + theme(plot.margin = unit(c(0,1,0,1), "cm"))
noleg <- theme(legend.position="none")
legsize <- theme(legend.title=element_text(size=10), legend.text=element_text(size=8))

plot.tsnes <- function(input){
  temp <- sum.data
  plotname <- input$color1
  
  if (input$color1 == "gene"){
     gene <- input$gene1
     temp <- cbind(sum.data,log2(R[gene,]+1))
     colnames(temp)<-c(colnames(sum.data),paste(gene, "\nlog2(RPKM+1)", sep=""))
     plotname <- paste(gene, "\nlog2(RPKM+1)", sep="")
  }

  p1 <- ggplot() + geom_jitter(data=temp,aes_(x=as.name("x.before"),y=as.name("y.before"),color=as.name(plotname)),na.rm=TRUE) + noax + noleg + ggtitle("Before grafting")
  p2 <- ggplot() + geom_jitter(data=temp,aes_(x=as.name("x.after"),y=as.name("y.after"),color=as.name(plotname)),na.rm=TRUE) + noax + legsize + ggtitle("After grafting")


  if (input$color1 == "gene" | plotname == "nGene"){
     by <- round(max(temp[,plotname])/6)
     if (plotname == "nGene") { by = 1000 }
     if (by < 1) {  by <- round(max(temp[,plotname])/0.6)/10 } 
     scale <- scale_color_gradientn(colors=c("green","yellow","red"), breaks =  seq(0, max(temp[,plotname])*1.1, by =
     by))
     p1 <- p1 + scale
     p2 <- p2 + scale	
  }else if (input$color1 == "Cluster") {
       scale <- scale_colour_manual(values = c("blue1", "red1","magenta", "green3", "cyan", "green","orange","yellow2"))
       p1 <- p1 + scale
       p2 <- p2 + scale
  }else if (input$color1 == "Source") {
       scale <- scale_colour_manual(values = c("grey56", "blue", "hotpink1", "red"))
       p1 <- p1 + scale
       p2 <- p2 + scale
  }else if (input$color1 == "SeuratCluster") {
       scale <- scale_colour_manual(values = c("green","orange","cyan","yellow2","green4", "deepskyblue","blue2","green4","green1","greenyellow","red1","magenta"))
       p1 <- p1 + scale
       p2 <- p2 + scale
  }				

  plot <- arrangeGrob(p1,p2,ncol=2,widths=c(3,4))
  class(plot) <- c("myplot", class(plot))
  return(plot)
}
  

######################################
# the actual server

shinyServer(  function(input, output,session) {

  plot <- reactive( { plot.tsnes(input) })
  output$plot_tsne <- renderPlot( { print(plot()) } )

  output$download <- downloadHandler(
    filename <- "tsnes_grafting.pdf",
    content = function(file) {
      pdf(file,width=10,height=4)              
      print(plot())
      dev.off()
  })  
})




