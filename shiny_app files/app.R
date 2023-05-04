library(shiny)
library(shinythemes)
library(ape)
library(BiocManager)
options(repos = BiocManager::repositories())
library(phytools)
library(devtools)
library(phangorn)
library(treeio)
library(paleotree)
library(geiger)
library(ggplot2)
library(phylobase)
library(data.table)
library(ggtree)
library(FossilSim)
library(shinyWidgets)
library(ggfortify)
library(shinyjs)
library(stats)
library(ggrepel)

ui <- fluidPage(theme=shinytheme("sandstone"),
  shinyjs::useShinyjs(),
  titlePanel(title=div(img(src="image.jpg"),"TaxonomR: Taxonomic Decision Support tool"), windowTitle = "TaxonomR" ), 
  sidebarLayout(
    sidebarPanel(
      fileInput("tre", "Upload tree file with node support", multiple=FALSE, accept=c(".tree", ".tre", ".txt")), 
      uiOutput("tre2"),
      div(style="display:inline-block;width:100%;text-align: right;",actionButton("add", "Add Tree")),
      fileInput("dat", "Upload data file", multiple=FALSE, accept=c("text/csv", "text/comma-separated-values, text/plain", ".csv")),
      #prettyCheckbox("combos", "Check possible group combinations", value = FALSE, shape = "curve", outline = FALSE, fill = FALSE, thick = FALSE, bigger = TRUE),
      sliderInput("mono", "Monophyly weight", value=0.0, min=0, max=1, step=0.05),
      sliderInput("sup", "Node support weight", value=0.0, min=0, max=1, step=0.05),
      sliderInput("br_ln", "Branch length weight", value=0.0, min=0, max=1, step=0.1), 
      sliderInput("val", "Current usage (stability) weight", value=0.0, min=0, max=1, step=0.05), 
      sliderInput("morph", "Morphological support weight", value=0.0, min=0, max=1, step=0.05),
      sliderInput("age", "Clade age weight", value=0.0, min=0, max=1, step=0.05),
      uiOutput("age_range"),
      sliderInput("geo", "Distribution weight", value=0.0, min=0, max=1, step=0.05),
      sliderInput("eco", "Ecological similarity weight", value=0.0, min=0, max=1, step=0.05),
      sliderInput("misc_1", "Custom trait/criterion 1", value=0.0, min=0, max=1, step=0.05),
      sliderInput("misc_2", "Custom trait/criterion 2", value=0.0, min=0, max=1, step=0.05)
   ), 
      mainPanel(
        tabsetPanel(id="tabs",
          tabPanel("Groups", 
                  plotOutput("tr"), 
                  fluidRow(
                  column(4, "", tableOutput("gen")), 
                  column(4, "", selectInput("hi_clade","Select clade to highlight", "")),
                  column(4, "", tableOutput("class"))
        )),
          tabPanel("Classifications", 
                   plotOutput("pcas"), 
                   fluidRow(h4("Please enter at least two groups to build a plot"),
                     column(3,wellPanel(
                       textInput("class_name", "Classification name", value="", placeholder = "enter unique label"),
                       uiOutput("checkGroup"),
                       actionButton("add_class", "Add to comparison"))), 
                     column(9,wellPanel(
                       tableOutput("classes_compare")
                     ))
        )
        
      )
          )
    )
  )
  
)
server <- function(input, output, session) {
 
 
############### ADD TREE ##########################
  observeEvent(input$add, {
    output$tre2 <- renderUI({
    fileInput("tre2", "Add tree with node ages", multiple=FALSE, accept=c(".tree", ".tre", ".txt"))
    })
  })
  
  observeEvent({
    input$tre
    input$dat
    },{ 
     output$age_range <- renderUI({
     sliderInput("age_range", "Age band range", min=0, max=lims()[1], value=c(lims()[2], lims()[3]))
    })
  })
  
  observeEvent({
    input$tre
    input$dat
  },{ 
    output$checkGroup <- renderUI({
      checkboxGroupInput("checkGroup", label = "Select groups", choices = dataset()$Group)
    })
  })

########################################################
#Diactivate slider input if the "Classification" tab is active and activate them back if the "Groups" tab is active
  observeEvent(
    input$tabs, {
      if(input$tabs=="Classifications") {
      shinyjs::disable("mono")
      shinyjs::disable("sup")
      shinyjs::disable("br_ln")
      shinyjs::disable("val")
      shinyjs::disable("morph")
      shinyjs::disable("age")
      shinyjs::disable("age_range")
      shinyjs::disable("geo")
      shinyjs::disable("eco")
	  shinyjs::disable("misc_1")
      shinyjs::disable("misc_2")
      }else{
        shinyjs::enable("mono")
        shinyjs::enable("sup")
        shinyjs::enable("br_ln")
        shinyjs::enable("val")
        shinyjs::enable("morph")
        shinyjs::enable("age")
        shinyjs::enable("age_range")
        shinyjs::enable("geo")
        shinyjs::enable("eco")
		shinyjs::enable("misc_1")
        shinyjs::enable("misc_2")
      }
     
  })
################################################
# CLASSIFICATION TAB INPUTS AND CALCULATIONS   #
################################################  
  # Read-in data file  
  dataset <- reactive({
    req(input$dat)
    d <- input$dat
    fread(d$datapath)
  })

  class_tab <- reactiveValues(df_class=NULL)
  
  ClassCompare <- eventReactive(input$add_class, {
    if(typeof(tree()[[2]])!="S4") {
      tree <- tree()[[2]]
    } else {
      tree <- tree()[[2]]@phylo
    }
    name <- as.character(dataset()$Group)
    genera <- list()
    for(i in 1:length(name)){
      tx <- as.character(gsub("\"", "", dataset()$Taxa[i]))
      tx <- as.character(gsub(" ", "", tx))
      taxa <- strsplit(tx, ",")[[1]]
      genera[[i]] <- taxa
    }
    names(genera) <- name
    uni <- genera
    n_gr <- length(input$checkGroup)
    selected <- as.vector(unlist(strsplit(input$checkGroup,",")),mode="list")
    
    ############# MONO METRIC ###################################################
    monos <- c()
    m=0
    for(i in 1:length(selected)) {
      index <- which(dataset()$Group == selected[[i]])
      if(is.monophyletic(tree, uni[[index]])==TRUE) {m=m+1} else {m=m+0}
    }
    mono_score <- m/length(selected)
    ############ SPLIT LUMP METRIC #############################################
    splits=0
    for(i in 1:length(selected)){
      index <- which(dataset()$Group == selected[[i]])
      y <- length(uni[[index]])
      splits=splits+y
    }
    split_score <- (splits/n_gr)/splits
    
    ########### GEN CAT METRIC ################################################
    #gen_cat=0
    #for(i in 1:length(selected)) {
    #  index <- which(dataset()$Group == selected[[i]])
    #  if(isTRUE(dataset()$Basic_category[[index]]=="Yes")) {gen_cat=gen_cat+1} else {gen_cat=gen_cat}
    #}
    #gen_cat_score <- gen_cat/length(selected)
    
    ########### BIOGEO METRIC #################################################
    biogeo <- c()
    for(i in 1:length(selected)) {
      index <- which(dataset()$Group == selected[[i]])
      biogeo <- c(biogeo, dataset()$Distribution[[index]])
    }
    biogeo_score <- length(unique(biogeo))/length(selected)
    
    ########## ECO METRIC #####################################################
    eco <- c()
    for(i in 1:length(selected)) {
      index <- which(dataset()$Group==selected[[i]])
      eco <- c(eco, dataset()$Ecology[[index]])
    }
    eco_score <- length(unique(eco))/length(selected)

    ########### Misc_1 METRIC #################################################
    #misc_1 <- c()
    #for(i in 1:length(selected)) {
    #  index <- which(dataset()$Group == selected[[i]])
    #  misc_1 <- c(misc_1, dataset()$Misc_1[[index]])
    #}
    #if(is.na(misc_1)==TRUE){misc_1_score=0} else{
    #misc_1_score <- length(unique(misc_1))/length(selected)
    #}
    ########### MISC_2 METRIC #####################################################
    #misc_2 <- c()
    #for(i in 1:length(selected)) {
    #  index <- which(dataset()$Group == selected[[i]])
    #  misc_2 <- c(misc_2, dataset()$Misc_2[[index]])
    #}
    #if(is.na(misc_2)==TRUE){misc_2_score=0} else{
    #  misc_2_score <- length(unique(misc_2))/length(selected)
    #}
    ########## AGE METRIC ####################################################
    age <- c()
    support <- c()
    for(i in 1:length(selected)) {
      index <- which(dataset()$Group==selected[[i]])
      if(length(uni[[index]])==1) {
        tp=which(tree$tip.label==uni[[index]])
        mrca = Ancestors(tree, tp, "parent")
      } else{
        mrca <- findMRCA(tree, tips=uni[[index]])
      }
      age_clade <- dateNodes(tree)[[mrca]]
      age <- c(age, age_clade)
      if (length(uni[[index]])==Ntip(tree)) {pp=1} else {
        if(typeof(tree()[[2]])!="S4") {
          if(is.monophyletic(tree()[[1]], uni[[index]])==TRUE) {
            pp <- as.numeric(tree()[[1]]$node.label[mrca-Ntip(tree)])
          } else {pp=0}
        } else {
          if(is.monophyletic(tree, uni[[index]])==TRUE) {
            n <- which(tree()[[2]]@data$node==mrca)
            pp <- tree()[[2]]@data$posterior[n]
          } else {pp=0}
          
        }
        if(pp>1){
          pp=pp/100
        }
      }
      support <- c(support, pp)
    }
    if(length(age)==1){
      age_score=0
    }else{
      age_score <- sd(age)
    }
    sup_score <- sum(support)/length(selected)
    
    
    #temp <- rbind(class_tab$df_class, data.frame("Classification"=input$class_name, "Group"=toString(input$checkGroup), "Monophyly"=mono_score, "Split/Lump"=split_score, "General category"=gen_cat_score, "Node support"=sup_score, "Age"=age_score, "Biogeo"=biogeo_score, "Eco"=eco_score))
    temp <- rbind(class_tab$df_class, data.frame("Classification"=input$class_name, "Group"=toString(input$checkGroup), "Monophyly"=mono_score, "Split/Lump"=split_score, "Node support"=sup_score, "Age"=age_score, "Biogeo"=biogeo_score, "Eco"=eco_score))
    
    class_tab$df_class <- temp
  })
  
  
  
# Read-in trees
tree <- reactive({
 req(input$tre)
  t <- input$tre
  if(is.null(input$tre2) == T){t2 <- input$tre} else {
    t2 <- input$tre2
  }
    
##################### FUNCTION TO READ DIFFERENT TREE FORMATS ##########################
  read.file <- function(file.name) {
      file <- try(read.beast(file.name), silent=TRUE)
      if(class(file)=="try-error") {
        file <- try(read.nexus(file.name), silent=TRUE)
        if (class(file)=="try-error") {
          file <- read.tree(file.name)
        }
      }
      file
    }
    
#################### READ-IN TREES #####################################################
    
    tree <- read.file(t$datapath)
    tree2 <- read.file(t2$datapath)
    trees <- list(tree, tree2)
    return(trees)
})
  
  lims <- reactive({
    req(input$dat)
    req(input$tre)
    if(typeof(tree()[[1]])!="S4") {
    max <- round(tree.max(tree()[[2]]), 2)
    lim <- 0.1*(max - dataset()$Age[1])
    } else {
    max <- round(tree.max(tree()[[2]]@phylo),2)
    lim <- 0.1*(tree.max(tree()[[2]]@phylo)-dataset()$Age[1])
    }
    lower <- dataset()$Age[1]-lim
    upper <- dataset()$Age[1]+lim
    lims <- c(max,lower, upper)
    return(lims)
    
  })
  
  combos <- reactive({
    req(input$dat)
    req(input$tre)
    data <- dataset()
    name <- as.character(data$Group)
    genera <- list()
    for(i in 1:length(name)){
      tx <- as.character(gsub("\"", "", data$Taxa[i]))
      tx <- as.character(gsub(" ", "", tx))
      taxa <- strsplit(tx, ",")[[1]]
      genera[[i]] <- taxa
    }
    names(genera) <- name
    ######################## DEFINING UNIQUE COMBINATIONS #####################
    # if(input$combos) {
    #   x <- c(names(genera))
    # 
    #   #find all possible combinations of specified groups
    #   combinations <- list()
    #   for(i in 1:length(x)) {
    #     a <- combn(x, i, simplify=F)
    #     combinations <- c(combinations, a)
    #   }
    # 
    #   #get taxic composition for each combination and remove duplicates
    # 
    #   unique_groups <- list()
    #   for(i in 1:length(combinations)) {
    #     gr <- c()
    #     for (j in 1:length(combinations[[i]])) {
    #       n <- which(names(genera)==combinations[[i]][j])
    #       gr<- c(gr, genera[[n]])
    #       gr <- unique(gr)
    #       gr <- c(gr)
    #     }
    #     unique_groups[[i]] <- gr
    #   }
    #   uno <- duplicated(unique_groups)
    #   unique_combo <- combinations[-c(which(uno==TRUE))]
    #   unique_combo_nums <- list()
    #   for(i in 1:length(unique_combo)){
    #     b <- c(genera[unique_combo[[i]]])
    #     d <- unique(unlist(b, use.names=F))
    #     unique_combo_nums[[i]] <- sort(d)
    #   }
    #   combo_names <- list()
    #   for(i in 1:length(unique_combo)){
    #     combo_names[i] <- toString(unique_combo[[i]])
    #   }
    #   names(unique_combo_nums) <- combo_names
    # 
    #   uno1 <- duplicated(unique_combo_nums)
    #   dups <- c()
    #   for(i in 1:length(uno1)){
    #     if(uno1[i]==TRUE) {dups <- c(dups, i)}
    #   }
    #   uni <- unique_combo_nums[-dups]
    #   opts <- c(names(uni))
    # }else{opts <- data$Group}
    opts <- data$Group
    return(opts)
  })

  
  
  weights <- reactive({
    req(input$age_range)
    req(input$dat)
    req(input$tre)
    m <- input$mono
    s <- input$sup
    b <- input$br_ln
    mr <- input$morph
    v <- input$val
    a <- input$age
    g <- input$geo
    e <- input$eco
	msc1 <- input$misc_1
	msc2 <- input$misc_2
    sum_score <- m+s+b+mr+v+a+g+e+msc1+msc2
    m1 <- m/sum_score
    s1 <- s/sum_score
    b1 <- b/sum_score
    mr1 <- mr/sum_score
    v1 <- v/sum_score
    a1 <- a/sum_score
    g1 <- g/sum_score
    e1 <- e/sum_score
	msc11 <- msc1/sum_score
	msc21 <- msc2/sum_score
    
    trees <- tree()
    data <- dataset()
######## Get group names and composition ###############
    name <- as.character(data$Group)
    genera <- list()
    for(i in 1:length(name)){
      tx <- as.character(gsub("\"", "", data$Taxa[i]))
      tx <- as.character(gsub(" ", "", tx))
      taxa <- strsplit(tx, ",")[[1]]
      genera[[i]] <- taxa
    }
    names(genera) <- name

####### Get Geo ranges ####################
    geos <- list()
    for(i in 1:length(name)){
      distr <- as.character(gsub("\"", "", data$Distribution[i]))
      distr <- as.character(gsub(" ", "", distr))
      dist <- strsplit(distr, ",")[[1]]
      geos[[i]] <- dist
    }
    names(geos) <- name
    
####### Get Eco information ###########
    ecos <- list()
    for(i in 1:length(name)){
      eco <- as.character(gsub("\"", "", data$Ecology[i]))
      eco <- as.character(gsub(" ", "", eco))
      ec <- strsplit(eco, ",")[[1]]
      ecos[[i]] <- ec
    }
    names(ecos) <- name
    
####### Get Misc_1 information ###########
    misc_tr1 <- list()
    for(i in 1:length(name)){
      msctr1 <- as.character(gsub(" ", "", data$Misc_1[i]))
      mctr1 <- strsplit(msctr1, ",")[[1]]
      misc_tr1[[i]] <- mctr1
    }
    names(misc_tr1) <- name

    
####### Get Misc_2 information ###########
    misc_tr2 <- list()
    for(i in 1:length(name)){
      msctr2 <- as.character(gsub(" ", "", data$Misc_2[i]))
      mctr2 <- strsplit(msctr2, ",")[[1]]
      misc_tr2[[i]] <- mctr2
    }
    names(misc_tr2) <- name
    
       
    
######################### DEFINING UNIQUE COMBINATIONS #####################
    # if(input$combos) {   
    #   x <- c(names(genera))
    #   
    #   #find all possible combinations of specified groups
    #   combinations <- list()
    #   for(i in 1:length(x)) {
    #     a <- combn(x, i, simplify=F)
    #     combinations <- c(combinations, a)
    #   }
    #   
    #   #get taxic composition for each combination and remove duplicates
    #   
    #   unique_groups <- list()
    #   for(i in 1:length(combinations)) {
    #     gr <- c()
    #     for (j in 1:length(combinations[[i]])) {
    #       n <- which(names(genera)==combinations[[i]][j])
    #       gr<- c(gr, genera[[n]])
    #       gr <- unique(gr)
    #       gr <- c(gr)
    #     }
    #     unique_groups[[i]] <- gr
    #   }
    #   
    #   
    #   uno <- duplicated(unique_groups)
    #   
    #   if(length(unique(uno))==1){
    #     unique_combo=combinations}else{
    #       unique_combo <- combinations[-c(which(uno==TRUE))]
    #     }
    #   
    #   
    #   unique_combo_nums <- list()
    #   for(i in 1:length(unique_combo)){
    #     b <- c(genera[unique_combo[[i]]])
    #     d <- unique(unlist(b, use.names=F))
    #     unique_combo_nums[[i]] <- sort(d)
    #   }
    #   
    #   combo_names <- list()
    #   for(i in 1:length(unique_combo)){
    #     combo_names[i] <- toString(unique_combo[[i]])
    #   }
    #   
    #   names(unique_combo_nums) <- combo_names
    #   
    #   uno1 <- duplicated(unique_combo_nums)
    #   dups <- c()
    #   for(i in 1:length(uno1)){
    #     if(uno1[i]==TRUE) {dups <- c(dups, i)}
    #   }
    #   uni <- unique_combo_nums[-dups]
    # } else {uni <- genera}
    uni <- genera
#################################################################################
    scores <- c()
    mrcas <- c()
###########################################################
if(typeof(trees[[1]])!="S4"){    
    
    for (i in 1:length(uni)) {
      
      score=0
      if(is.monophyletic(trees[[1]], uni[[i]])==TRUE) {score=score+m1} else {score=score}
  ###################### ADD CHECK FOR IF THE GROUP IS MONOTYPIC#######################
      if(length(uni[[i]])==1) {
        tp=which(trees[[1]]$tip.label==uni[[i]])
        mrca = Ancestors(trees[[1]], tp, "parent")
        } else{
        mrca <- findMRCA(trees[[1]], tips=uni[[i]])
        }
      mrcas <- c(mrcas, mrca)
      if(is.monophyletic(trees[[1]], uni[[i]])==TRUE) {
      pp <- as.numeric(trees[[1]]$node.label[mrca-Ntip(trees[[1]])])
      score = score+0.01*pp*s1
      } else {score=score}
      
      
      br <- which.edge(trees[[2]], uni[[i]])
      if(length(br)==1 || length(br)==2){
        ratio <- max(trees[[2]]$edge.length[br])/trees[[2]]$edge.length[br[1]-1]
      } else if(length(br)==Nedge(trees[[2]])){
        #ratio <- max(trees[[2]]$edge.length[br])/trees[[2]]$root.edge #any number greater than 1 to result in 0 additional score
        ratio <- 0.5
      }
      else {
        br1 <- br[ which(trees[[2]]$edge[br,][,2]>Ntip(trees[[2]]))]
        ratio <- max(trees[[2]]$edge.length[br1])/trees[[2]]$edge.length[br[1]-1]
        if(length(ratio)==0) {
          ratio <- max(trees[[2]]$edge.length[br])/trees[[2]]$root.edge
        }
      }
      
      if(ratio > 1) {score = score} else if(ratio == 1) {score=score+0.5*b1} else if(ratio < 1) {score = score+b1}
      
      usg <- which(data$Group==names(uni[i]))
      if(isTRUE(data$Usage[usg]=="major_use")) {score=score+v1} else if(isTRUE(data$Usage[usg]=="minor_use"|| data$Usage[usg]=="new_use")) {score=score+0.5*v1} else if(isTRUE(data$Usage[usg]=="new")) {score=score} else {score=score} 
      
      
      if(isTRUE(data$Morph[usg]=="synap")) {score=score+mr1} else if(isTRUE(data$Morph[usg]=="combo")) {score=score+0.5*mr1} else if(isTRUE(data$Morph[usg]=="plastic")) {score=score+0.25*mr1} else if(isTRUE(data$Morph[usg]=="none")) {score=score} else {score=score}
      
      date <- dateNodes(trees[[2]])[[mrca]]
      if(date > input$age_range[1] && date < input$age_range[2]) {score=score+a1} else {score=score}
      
      if(length(geos[[i]])==1){score= score+g1} else {score=score}
      
      if(length(ecos[[i]])==1){score=score+e1} else{score=score}
      if(length(misc_tr1[[i]])==1){score=score+msc11} else{score=score}
      if(length(misc_tr2[[i]])==1){score=score+msc21} else{score=score}

      scores <- c(scores, score)
      
    }
} else {
  for (i in 1:length(uni)) {
    
    score=0
    if(is.monophyletic(trees[[1]]@phylo, uni[[i]])==TRUE) {score=score+m1} else {score=score}
    
    if(length(uni[[i]])==1) {
      tp=which(trees[[1]]@phylo$tip.label==uni[[i]])
      mrca = Ancestors(trees[[1]]@phylo, tp, "parent")
    } else{
      mrca <- findMRCA(trees[[1]]@phylo, tips=uni[[i]])
    }
    mrcas <- c(mrcas, mrca)
    
    if(is.monophyletic(trees[[1]]@phylo, uni[[i]])==TRUE) {
      n <- which(trees[[1]]@data$node==mrca)
      pp <- trees[[1]]@data$posterior[n]
      if(pp < 1){score = score+pp*s1} else{
      score = score+0.01*pp*s1}
    } else {score=score}
    

    
    br <- which.edge(trees[[1]]@phylo, uni[[i]])
    if(length(br)==1 || length(br)==2){
      ratio <- max(trees[[1]]@phylo$edge.length[br])/trees[[1]]@phylo$edge.length[br[1]-1]
    } else if(length(br)==Nedge(trees[[1]]@phylo)){
      #ratio <- max(trees[[1]]@phylo$edge.length[br])/trees[[1]]@phylo$root.edge #any number greater than 1 to result in 0 additional score
      ratio <- 0.5
    }
    else {
      br1 <- br[which(trees[[1]]@phylo$edge[br,][,2]>Ntip(trees[[1]]@phylo))]
      ratio <- max(trees[[1]]@phylo$edge.length[br1])/trees[[1]]@phylo$edge.length[br[1]-1]
    }
    
   if(ratio > 1) {score = score} else if(ratio == 1) {score=score+0.5*b1} else if(ratio < 1) {score = score+b1}
    
    
    usg <- which(data$Group==names(uni[i]))
    if(isTRUE(data$Usage[usg]=="major_use")) {score=score+v1} else if(isTRUE(data$Usage[usg]=="minor_use"|| data$Usage[usg]=="new_use")) {score=score+0.5*v1} else if(isTRUE(data$Usage[usg]=="new")) {score=score} else {score=score} 
    
    
    if(isTRUE(data$Morph[usg]=="synap")) {score=score+mr1} else if(isTRUE(data$Morph[usg]=="combo")) {score=score+0.5*mr1} else if(isTRUE(data$Morph[usg]=="plastic")) {score=score+0.25*mr1} else if(isTRUE(data$Morph[usg]=="none")) {score=score} else {score=score}
    
    date <- dateNodes(trees[[1]]@phylo)[[mrca]]
    if(date > input$age_range[1] && date < input$age_range[2]) {score=score+a1} else {score=score}
    
    if(length(geos[[i]])==1){score= score+g1} else {score=score}
    
    if(length(ecos[[i]])==1){score=score+e1} else{score=score}
    if(length(misc_tr1[[i]])==1){score=score+msc11} else{score=score}
    if(length(misc_tr2[[i]])==1){score=score+msc21} else{score=score}

    
    scores <- c(scores, score)
    
  }
}
    score_table <- data.frame(group=names(uni), score_out_of_1.0=scores)
    score_table <- score_table[order(-scores),]
    score_table
    
  })
  

  observe({
    updateSelectInput(session, "hi_clade",
                      label = "Select Group",
                      choices=combos(),
                      selected=NULL
    ) 
    
  })
  
 hilight <- reactive({
  req(input$tre)
  #req(input$tre2)
  req(input$dat)
  req(input$hi_clade)
  tree1 <- tree()  
  data1 <- dataset()
  data <- dataset()
  
  name <- as.character(data1$Group)
  genera <- list()
  for(i in 1:length(name)){
    tx <- as.character(gsub("\"", "", data$Taxa[i]))
    tx <- as.character(gsub(" ", "", tx))
    taxa <- strsplit(tx, ",")[[1]]
    genera[[i]] <- taxa
  }
  names(genera) <- name
  ######################## DEFINING UNIQUE COMBINATIONS #####################
  # if(input$combos) { 
  #   x <- c(names(genera))
  #   
  #   #find all possible combinations of specified groups
  #   combinations <- list()
  #   for(i in 1:length(x)) {
  #     a <- combn(x, i, simplify=F)
  #     combinations <- c(combinations, a)
  #   }
  #   
  #   #get taxic composition for each combination and remove duplicates
  #   
  #   unique_groups <- list()
  #   for(i in 1:length(combinations)) {
  #     gr <- c()
  #     for (j in 1:length(combinations[[i]])) {
  #       n <- which(names(genera)==combinations[[i]][j])
  #       gr<- c(gr, genera[[n]])
  #       gr <- unique(gr)
  #       gr <- c(gr)
  #     }
  #     unique_groups[[i]] <- gr
  #   }
  #   uno <- duplicated(unique_groups)
  #   unique_combo <- combinations[-c(which(uno==TRUE))]
  #   unique_combo_nums <- list()
  #   for(i in 1:length(unique_combo)){
  #     b <- c(genera[unique_combo[[i]]])
  #     d <- unique(unlist(b, use.names=F))
  #     unique_combo_nums[[i]] <- sort(d)
  #   }
  #   combo_names <- list()
  #   for(i in 1:length(unique_combo)){
  #     combo_names[i] <- toString(unique_combo[[i]])
  #   }
  #   names(unique_combo_nums) <- combo_names
  #   
  #   uno1 <- duplicated(unique_combo_nums)
  #   dups <- c()
  #   for(i in 1:length(uno1)){
  #     if(uno1[i]==TRUE) {dups <- c(dups, i)}
  #   }
  #   uni <- unique_combo_nums[-dups]
  # } else {uni <- genera}
  uni <- genera
  
  k <- which(names(uni)==input$hi_clade) 
  #tx <- as.character(gsub("\"", "", dataset()$Taxa[k]))
  #tx <- as.character(gsub(" ", "", tx))
  #taxa <- strsplit(tx, ",")[[1]]
  taxa <- uni[[k]]
  
if(typeof(tree()[[2]])!="S4"){
  if(length(taxa)==1){
    e <- rep("black", Nedge(tree()[[2]]))
    e[which.edge(tree()[[2]], taxa)] <- "red"
    tc <- rep("black", Ntip(tree()[[2]]))
    tc[which(tree()[[2]]$tip.label==taxa)] <- "red"
    plot(tree()[[2]], edge.color=e, tip.color=tc, no.margin=T, font=1, cex=0.5, edge.width=1.5)
  } else {
    zoom(tree()[[2]], taxa, subtree = TRUE, no.margin = TRUE)
  }
} else {
  if(length(taxa)==1){
    e <- rep("black", Nedge(tree()[[2]]@phylo))
    e[which.edge(tree()[[2]]@phylo, taxa)] <- "red"
    tc <- rep("black", Ntip(tree()[[2]]@phylo))
    tc[which(tree()[[2]]@phylo$tip.label==taxa)] <- "red"
    plot(tree()[[2]]@phylo, edge.color=e, tip.color=tc, no.margin=T, font=1, cex=0.5, edge.width=1.5)
  } else {
  zoom(tree()[[2]]@phylo, taxa, subtree = TRUE, no.margin = TRUE)
  }
}
  })
 
classification <- reactive({
  req(input$age_range)
  req(input$dat)
  req(input$tre)
  req(weights())
  trees <- tree()
  data <- dataset()
  ######## Get group names ###############
  name <- as.character(data$Group)
  genera <- list()
  for(i in 1:length(name)){
    tx <- as.character(gsub("\"", "", data$Taxa[i]))
    tx <- as.character(gsub(" ", "", tx))
    taxa <- strsplit(tx, ",")[[1]]
    genera[[i]] <- taxa
  }
  names(genera) <- name
  # if(input$combos)  {
  #   x <- c(names(genera))
  #   combinations <- list()
  #   for(i in 1:length(x)) {
  #     a <- combn(x, i, simplify=F)
  #     combinations <- c(combinations, a)}
  #   unique_groups <- list()
  #   for(i in 1:length(combinations)) {
  #     gr <- c()
  #     for (j in 1:length(combinations[[i]])) {
  #       n <- which(names(genera)==combinations[[i]][j])
  #       gr<- c(gr, genera[[n]])
  #       gr <- unique(gr)
  #       gr <- c(gr)
  #     }
  #     unique_groups[[i]] <- gr
  #   }
  #   uno <- duplicated(unique_groups)
  #   unique_combo <- combinations[-c(which(uno==TRUE))]
  #   unique_combo_nums <- list()
  #   for(i in 1:length(unique_combo)){
  #     b <- c(genera[unique_combo[[i]]])
  #     d <- unique(unlist(b, use.names=F))
  #     unique_combo_nums[[i]] <- sort(d)
  #   }
  #   
  #   combo_names <- list()
  #   for(i in 1:length(unique_combo)){
  #     combo_names[i] <- toString(unique_combo[[i]])
  #   }
  #   
  #   names(unique_combo_nums) <- combo_names
  #   
  #   uno1 <- duplicated(unique_combo_nums)
  #   dups <- c()
  #   for(i in 1:length(uno1)){
  #     if(uno1[i]==TRUE) {dups <- c(dups, i)}
  #   }
  #   uni <- unique_combo_nums[-dups]
  # } else {uni <- genera}
  uni <- genera
  if(typeof(tree()[[2]])!="S4") {
    taxa <- tree()[[2]]$tip.label
  }else {
    taxa <- tree()[[2]]@phylo$tip.label
  }
  
  suggested_df <- data.frame(Names=character(), Score=double())
  for(x in 1:length(uni)){
    if(typeof(tree()[[2]])!="S4") {
      taxa <- tree()[[2]]$tip.label
    }else {
      taxa <- tree()[[2]]@phylo$tip.label
    }
    classes <- c()
    points <- c()
    taken_x <- uni[[as.numeric(row.names(weights())[x])]]
    pts_x <- weights()[x,2]
    taxa <- setdiff(taxa, taken_x)
    for(j in 1:length(uni)) {	
      if(j == x && j < length(uni)) {j=j+1}
      taken <- uni[[as.numeric(row.names(weights())[j])]]
      pts <- weights()[j,2]
      overlap <- setdiff(taken, taxa)
      if(length(overlap)>0) {print ('skipping this group'); next}
      left <- setdiff(taxa, taken)
      taxa <- left
      classes <- c(classes, j)
      points <- c(points, pts)
      
    }
    classes <- c(classes, x)
    points <- c(points, pts_x)
    fin_pts = sum(points)/length(points)
    classification_x <- c()
    for(y in 1:length(classes)) {
      classification_x <- c(classification_x, names(uni[as.numeric(row.names(weights())[classes[y]])]))
      classification_x <- sort(classification_x)
    }
    df = data.frame(Names=toString(classification_x), Score=fin_pts)
    suggested_df = rbind(suggested_df, df)
  } 


  suggested_df1 <- data.frame(Names=character(), Score=double())
  for(x in length(uni):1){
    if(typeof(tree()[[2]])!="S4") {
      taxa <- tree()[[2]]$tip.label
    }else {
      taxa <- tree()[[2]]@phylo$tip.label
    }
    classes <- c()
    points <- c()
    taken_x <- uni[[as.numeric(row.names(weights())[x])]]
    pts_x <- weights()[x,2]
    taxa <- setdiff(taxa, taken_x)
    for(j in length(uni):1) {	
      if(j == x && j < length(uni)) {j=j+1}
      taken <- uni[[as.numeric(row.names(weights())[j])]]
      pts <- weights()[j,2]
      overlap <- setdiff(taken, taxa)
      if(length(overlap)>0) {print ('skipping this group'); next}
      left <- setdiff(taxa, taken)
      taxa <- left
      classes <- c(classes, j)
      points <- c(points, pts)
      
    }
    classes <- c(classes, x)
    points <- c(points, pts_x)
    fin_pts = sum(points)/length(points)
    classification_x <- c()
    for(y in 1:length(classes)) {
      classification_x <- c(classification_x, names(uni[as.numeric(row.names(weights())[classes[y]])]))
      classification_x <- sort(classification_x)
    }
    df = data.frame(Names=toString(classification_x), Score=fin_pts)
    suggested_df1 = rbind(suggested_df1, df)
  } 
  
  suggested_df1a <- data.frame(Names=character(), Score=double())
  mid_count <- seq(from=length(uni)-3, to=1)
  print(mid_count)
  mid_count2 <- c(length(uni):(length(uni)-2))
  mid_count <- c(mid_count, mid_count2)
  print(mid_count)
  for(x in mid_count){
    if(typeof(tree()[[2]])!="S4") {
      taxa <- tree()[[2]]$tip.label
    }else {
      taxa <- tree()[[2]]@phylo$tip.label
    }
    md=length(uni) - 2
    classes <- c()
    points <- c()
    taken_x <- uni[[as.numeric(row.names(weights())[x])]]
    pts_x <- weights()[x,2]
    taxa <- setdiff(taxa, taken_x)
    for(j in mid_count) {	
      if(j == x && j < length(uni)) {j=j+1}
      taken <- uni[[as.numeric(row.names(weights())[j])]]
      pts <- weights()[j,2]
      overlap <- setdiff(taken, taxa)
      if(length(overlap)>0) {print ('skipping this group'); next}
      left <- setdiff(taxa, taken)
      taxa <- left
      classes <- c(classes, j)
      points <- c(points, pts)
      
    }
    classes <- c(classes, x)
    points <- c(points, pts_x)
    fin_pts = sum(points)/length(points)
    classification_x <- c()
    for(y in 1:length(classes)) {
      classification_x <- c(classification_x, names(uni[as.numeric(row.names(weights())[classes[y]])]))
      classification_x <- sort(classification_x)
    }
    df = data.frame(Names=toString(classification_x), Score=fin_pts)
    suggested_df1a = rbind(suggested_df1a, df)
  } 
  
  suggested_df2 = rbind(suggested_df, suggested_df1, suggested_df1a)
    
  classification <- suggested_df2[order(suggested_df2$Score,decreasing=TRUE),]
  classification_out <- classification[!duplicated(classification), ]
  #classification_out <- classification
  
  # classification <- c()
  #for(y in 1:length(classes)) {
  #  classification <- c(classification, names(uni[as.numeric(row.names(weights())[classes[y]])]))
  #}
  #classification <- c(classification, fin_pts)
  
  # classification_out <- data.frame("Suggested classification" = classification)
  classification_out
}) 


dataPCA <- reactive({
  req(input$tre)
  req(input$dat)
  req(ClassCompare())
  classification_df <- as.data.frame(ClassCompare())
  type <- sub(".*_", "", classification_df$Classification)
  classification_df$Type <- type
  if(nrow(classification_df)>=2) {
    active_data <- classification_df[,3:8]
    row.names(active_data) <- classification_df[,1]
    class_pca <- prcomp(active_data, center = T, scale. = T)
    pca_loads <- as.data.frame(class_pca$rotation)
    return(pca_loads)
  }
})

plotPCA <- reactive({
  req(input$tre)
  req(input$dat)
  req(ClassCompare())
  classification_df <- as.data.frame(ClassCompare())
  type <- sub(".*_", "", classification_df$Classification)
  classification_df$Type <- type
  if(nrow(classification_df)>=2) {
    active_data <- classification_df[,3:8]
    row.names(active_data) <- classification_df[,1]
    autoplot(prcomp(active_data), data=classification_df, colour="Type", label=F, size=7) + geom_text_repel(label=row.names(active_data), size=6)+theme(legend.text=element_text(size=16))
    
  }
})

  
  output$tr <- renderPlot({
    hilight()
  }) 

  
  output$gen <- renderTable({
    weights()
  })
  
  output$class <- renderTable({
    classification()  
  })
  
  output$pcas <- renderPlot({
    plotPCA()
  })
  
  output$classes_compare <- renderTable({
    ClassCompare()
  })
}

shinyApp(ui = ui, server = server)