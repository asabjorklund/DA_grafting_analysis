## Written 2017 by Asa Bjorklund asa.bjorklund@scilifelab.se

options <- list("Cell source"="Source","Cluster"="Cluster","Cluster Seurat"="SeuratCluster","Detected Genes"= "nGene", "Gene Expression"="gene")
Gene_names <- rownames(R)
genes <- as.data.table(Gene_names)


# Layout UI custom components:
selectColorBy_wrapper <- function(id, sel) {
  selectInput(id, "Color cells by (also affects boxplots below):",
              choices=colorChoices, selected=sel)
}
selectGene_wrapper <- function(id,sel='') {
  selectInput(paste("gene",id,sep=""), "Select gene:",
                 choices=c(Start_typing = '', genes),
                 selected=sel)
}



shinyUI(fluidPage(
  title = "Grafting data SS2",
  wellPanel(style="background-color: #ffffff;",
  h1("Grafting data SS2"),
  p("Visualization of data from single cell RNA sequencing with Smartseq2 protocol of fetal and clinical grade ventral midbrain (VM)-patterned human embryonic stem cells (hESCs) before and after long-term survival and functional maturation in a pre-clinical rat model of PD."),
  p("Please select what type of data you would like to show for each of the two t-SNE layouts below (before and after transplantation). When selecting \"Gene expression\" you also need to specify the gene in the lower box. To show your selections, click on \"Apply Changes\". Plots can also be exported as pdf with the \"Save as pdf\" button.")
  ),
  
  wellPanel(style="background-color: #ffffff;",

  fluidRow(
      column(12,
        plotOutput("plot_tsne")
      )
  ), 

  fluidRow(
    column(3,offset=3,
      h2("Select colouring"),
      selectInput("color1", label = h3("Color cells by"),
              choices = options, selected = "Cluster"),
      selectGene_wrapper(1),	      
      helpText("Note: \"Color cells by\" must be set to \"Gene Expression\" to color by expression")	
    ),
    column(2,offset=2,
	br(),
        submitButton("Apply Changes"),       
	br(),
        downloadButton('download', 'Save as pdf')
    )
  ),
  br(),
  p("")
)))

