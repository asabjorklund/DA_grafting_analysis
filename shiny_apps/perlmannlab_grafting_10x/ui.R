## Written 2017 by Asa Bjorklund asa.bjorklund@scilifelab.se

options <- list("Cluster"="CelltypeCluster","Sample"="orig.ident","Detected Genes"= "nFeature_RNA", "Gene Expression"="gene", "FACS sorting"="Type","Rat"="Rat")


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
  title = "Grafting data 10x",
  wellPanel(style="background-color: #ffffff;",
  h1("Grafting data 10x"),
  p("Visualization of data from single cell RNA sequencing using 10x Genomics protocol of clinical grade ventral midbrain (VM)-patterned human embryonic stem cells after long-term survival and functional maturation in a pre-clinical rat model of PD."),
  p("Please select what type of data you would like to show for each of the UMAP or t-SNE layout. When selecting \"Gene expression\" you also need to specify the gene in the lower box. To show your selections, click on \"Apply Changes\". Plots can also be exported as pdf with the \"Save as pdf\" button.")
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
      selectInput("color1", label = h3("Color cells by:"),
              choices = options, selected = "Cluster"),
      selectGene_wrapper(1),	      
      helpText("Note: \"Color cells by\" must be set to \"Gene Expression\" to color by expression")	
    ),
    column(2,offset=2,
	br(),
        submitButton("Apply Changes"),       
	br(),
        downloadButton('download', 'Save as pdf'),
	br(),
	selectInput("reduction", label = h3("Dimensionality reduction"),
              choices = list("UMAP"="UMAP","tSNE"="tSNE"), selected = "UMAP")
    )
  ),
  br(),
  p("")
)))

