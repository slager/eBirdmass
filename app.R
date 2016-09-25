library(shiny)
library(rvest)
library(httr)
library(DT)

## Overarching function to scrape eB checklist and calculate outputs
maketable <- function(URL){
  ## Grab S# from user input and standardize URL
  regmatches(URL,regexpr("S[0123456789]+",URL,perl=TRUE)) -> S
  paste0("http://ebird.org/ebird/view/checklist/",S) -> URL

  ## Scrape eBird checklist
  #xhtml <- read_html(URL) #Simplest but throws error:
  #"Peer certificate cannot be authenticated with given CA certificates"
  #http://stackoverflow.com/questions/34551299/how-to-pass-ssl-verifypeer-in-rvest
  set_config(config(ssl_verifypeer = 0L)) #Don't verify certificate
  xhtml <- read_html(content(GET(URL), as="text")) #httr workaround for read_html
  species<-html_nodes(x=xhtml,css=".se-name") %>% html_text()
  count<-html_nodes(x=xhtml,css=".se-count") %>% html_text()
  title<-html_nodes(x=xhtml,css="title") %>% html_text()
  title<-substr(title,19,nchar(title)) #remove "eBird Checklist - " from title
  
  #Evaluate warning for X counts and coerce to 1s
  ifX <- FALSE
  if ("X" %in% count){TRUE->ifX ; 1 -> count[which(count=="X")]}
  as.numeric(count) -> count
  
  #Import eBird taxonomy and masses files
  tax <- read.csv("eBird_Taxonomy_v2016_9Aug2016.csv",header=T,stringsAsFactors=F,na.strings="")
  m <- read.csv("masses_tax2016.csv",header=T,stringsAsFactors=F,na.strings=c("","NA"))

  #Create working data frame d
  d <- data.frame("com1"=species,"com2"=NA,"count"=count,stringsAsFactors=F)
  
  #Add CATEGORY to main data frame, from Taxonomy file
  sapply(species,function(x){tax$CATEGORY[which(tax$PRIMARY_COM_NAME==x)]},USE.NAMES=F) -> d$CATEGORY
  
  #Add com2 (=taxon for mass matching) to d for each com1 (=original taxon)
  #Do this differently depending on CATEGORY
  apply(d[,c('com1','CATEGORY')],1,function(x){
    if(x['CATEGORY'] %in% c('species','spuh')){return(x['com1'])}
    if(x['CATEGORY'] %in% c('domestic','issf','form','intergrade')){
      which(strsplit(x['com1'],'')[[1]]=='(') -> p
      if (length(p)>0){return(substr(x['com1'],1,p-2))}}
    if(x['CATEGORY']=='hybrid'){
      match(x['com1'],tax$PRIMARY_COM_NAME) -> i
      regexpr(" [(]", tax$SCI_NAME[i])[[1]][1]-1 -> p
      regexpr(" x ", tax$SCI_NAME[i])[[1]][1]-1 -> q
      if(max(p,q)>0){substr(tax$SCI_NAME[i],1,max(p,q)) -> s}
      match(s,tax$SCI_NAME) -> i
      return(tax$PRIMARY_COM_NAME[i])}
    if(x['CATEGORY']=='slash'){
      match(x['com1'],tax$PRIMARY_COM_NAME) -> i
      regexpr("/", tax$SCI_NAME[i])[[1]][1]-1 -> p
      if(length(p)>0){substr(tax$SCI_NAME[i],1,p) -> s}
      match(s,tax$SCI_NAME) -> i
      return(tax$PRIMARY_COM_NAME[i])}
  }) -> d$com2
  
  #Add GET_MASS_CODE for com2 to d (later use for matching masses)
  sapply(d$com2,function(x){tax$SPECIES_CODE[which(tax$PRIMARY_COM_NAME==x)]},USE.NAMES=F) -> d$GET_MASS_CODE
  
  #Add mass & mass comments to d, using SPECIES_CODE
  sapply(d$GET_MASS_CODE,function(x){
    ifelse(x %in% m$SPECIES_CODE,m[which(m$SPECIES_CODE==x),'mass'],NA)
    }) -> d$mass
  sapply(d$GET_MASS_CODE,function(x){
    ifelse(x %in% m$SPECIES_CODE,m[which(m$SPECIES_CODE==x),'masscomm'],NA)
    }) -> d$masscomm
  "" -> d$masscomm[which(is.na(d$masscomm))]
  
  #Count taxa with NA masses
  sum(is.na(d$mass)) -> countNAmass
  
  #Calculate mass stats and summary stats
  d$count * d$mass -> d$subtotal #mass subtotal by taxon
  d$subtotal/sum(na.omit(d$subtotal)) -> d$pm # % total mass by taxon
  d$count/sum(na.omit(d$count)) -> d$pc # % total count by taxon

  ## Make function output data frame f
  f <- data.frame(stringsAsFactors=F,i=as.numeric(rownames(d)),taxon=d$com1,count=d$count,count.p=d$pc,unit.g=d$mass,subtotal.g=d$subtotal,mass.p=d$pm,note=d$masscomm)
  
  ## Standalone calculations to return from function
  sum(d$count) -> ni # Total count of all individuals
  sum(na.omit(d$subtotal)) -> tm # Total mass (g) of all individuals
  length(d$com1) -> ns # Total number of taxa
  
  return(list('countNAmass'=countNAmass,'title'=title,'URL'=URL,'tm'=tm,'ni'=ni,'ns'=ns,'ifX'=ifX,'f'=f))
}  ## End of giant maketable() function

## Begin Shiny UI
ui <- fluidPage(
  titlePanel("eBird Biomass Calculator"),
  fluidRow(
    column(3, wellPanel(
      textInput("text", "eBird Checklist URL:", ""),
      actionButton("do", "Calculate"),
      tags$br(),
      tags$br(),
      p(tags$a(href="https://github.com/slager/eBirdmass", "R code", target="_blank")," by ",tags$a(href="https://twitter.com/dlslager", "Dave Slager", target="_blank")),
      p("Masses",tags$a(href="https://www.amazon.com/Sibley-Guide-Birds-David-Allen/dp/0679451226","(Sibley 2000)", target="_blank"),"transcribed by Sean Fitzgerald")
    )),
    column(6,
           h5(htmlOutput("mySite")),
           DT::dataTableOutput("table"),
           h5(htmlOutput("coerceX")),
           h5(htmlOutput("massNA")),
           h5(htmlOutput("nins")),
           h5(htmlOutput("metric")),
           h5(htmlOutput("english")),
           h5(htmlOutput("pennies"))
    ) #End Column
  ) #End Fluidrow
) #End Fluidpage

## Begin Shiny Server
server <- function(input, output){
  eventReactive(input$do,{maketable(input$text)}) -> maketablereactive
  output$table <- DT::renderDataTable(DT::datatable({maketablereactive()[['f']]},
                                                    rownames=FALSE,
                                                    colnames=c("","Taxon","Count","%Count","Unit(g)","Subtotal(g)","%Mass","Note"),
                                                    options=list(paging=FALSE,searching=FALSE,info=FALSE)
                                                    ) %>%
                                        formatCurrency(columns=c(3,6),currency="",interval=3,mark=",",digits=0) %>%
                                        formatCurrency(columns=5,currency="",interval=3,mark=",",digits=1) %>%
                                        formatPercentage(c(4,7), 2))
  output$mySite <- renderUI({
    tags$a(href = maketablereactive()[['URL']], maketablereactive()[['title']], target="_blank")})
  output$massNA <- renderUI({
    tags$p(style="color:red",ifelse(maketablereactive()[['countNAmass']]>0,paste0("Warning: Mass=NA for ",maketablereactive()[['countNAmass']]," taxa"),""))})
  output$coerceX <- renderUI({
    tags$p(style="color:red",ifelse(maketablereactive()[['ifX']],"Warning: Xs were coerced to 1s",""))})
  output$nins <- renderUI({
    tags$p(paste0(format(maketablereactive()[['ni']],big.mark=",")," individuals"),
           tags$br(),
           paste0(maketablereactive()[['ns']]," taxa"))})
  output$metric <- renderUI({
    tags$p(paste0(format(round(maketablereactive()[['tm']],1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," grams"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]/1000,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," kg"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]/1e6,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," metric tonnes"))})
  output$english <- renderUI({
    tags$p(paste0(format(round(maketablereactive()[['tm']]*0.035274,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," oz"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]*0.00220462,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," lb"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]*1.10231e-6,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," tons"))})
  output$pennies <- renderUI({
    tags$p(paste0("$",format(round(maketablereactive()[['tm']]/2.5*.01, 2), nsmall = 2, big.mark=",")," in pennies"))})
} # End of Shiny Server
shinyApp(ui, server)