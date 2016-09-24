library(shiny)
library(rvest)
library(httr)
library(DT)

## Define function to scrape, calculate, and generate output data
maketable <- function(x){
  x -> URL
  ## Standardize URL for either of the 2 types or just S#
  if (grepl("checklist?subID=S",URL)){
    regexpr(pattern ='subID=S',URL) -> i
    substr(URL,i,nchar(URL)) -> S
    paste("http://ebird.org/ebird/view/checklist?",S,sep="") -> URL}
  if (grepl("checklist/S",URL)){
    regexpr(pattern ='/S',URL) -> i
    substr(URL,i+1,nchar(URL)) -> S
    paste("http://ebird.org/ebird/view/checklist?subID=",S,sep="") -> URL}
  if (substr(URL,1,1)=="S"){
    paste("http://ebird.org/ebird/view/checklist?subID=",URL,sep="") -> URL}
    
  ## Scrape eBird Checklist - numbered HTML lines
  #This method no longer worked when I updated R and rvest
  #checklist <- paste(capture.output(html(URL), file=NULL))

  #xhtml <- read_html(URL)
  # Causes error "Peer certificate cannot be authenticated with given CA certificates"
  #Workaround:
  #http://stackoverflow.com/questions/34551299/how-to-pass-ssl-verifypeer-in-rvest
  set_config(config(ssl_verifypeer = 0L))
  xhtml <- read_html(content(GET(URL), as="text"))
  
  species<-html_nodes(x=xhtml,css=".se-name") %>% html_text()
  count<-html_nodes(x=xhtml,css=".se-count") %>% html_text()
  title<-html_nodes(x=xhtml,css="title") %>% html_text()
  title<-substr(title,19,nchar(title))
  
  ##Grep and trim the species (This longer needed using html_nodes method)
  # grep('class="se-name"',checklist) -> i #grep lines of the common names
  # checklist[i] -> species
  # gsub("\t","",species) -> species #trim tabs
  # sub('<h5 class=\"se-name\">',"",species) -> species #trim left tags
  # sub('</h5>',"",species) -> species #trim right
  # #Grep and trim the counts
  # grep('class="se-count"',checklist) -> j #grep lines of the counts
  # checklist[j] -> count
  # sub('<th><h5 class=\"se-count\">',"",count) -> count #trim left
  # sub('</h5></th>',"",count) -> count #trim right
  #Evaluate warning for Xs and set Xs = 1
  FALSE -> ifX
  if ("X" %in% count){
    TRUE->ifX
    1 -> count[which(count=="X")]
  }
  #Convert counts to values
  as.numeric(count) -> count
  
  #Create data frame
  data.frame("com1"=species,"com2"=NA,"count"=count,stringsAsFactors=F) -> d
  
  #Get the eBird taxonomy
  read.csv("eBird_Taxonomy_v2016_9Aug2016.csv",header=T,stringsAsFactors=F,na.strings="") -> tax
  
  #Add CATEGORY to main data frame
  sapply(species,function(x){tax$CATEGORY[which(tax$PRIMARY_COM_NAME==x)]},USE.NAMES=F) -> d$CATEGORY
  
  #Add com2 to data frame for different levels of CATEGORY
  apply(d[,c('com1','CATEGORY')],1,function(x){
    if(x['CATEGORY'] %in% c('species','spuh')){return(x['com1'])}
    if(x['CATEGORY'] %in% c('domestic','issf','form','intergrade')){
      which(strsplit(x['com1'],'')[[1]]=='(') -> p
      if (length(p)>0){return(substr(x['com1'],1,p-2))}}
    if(x['CATEGORY']=='hybrid'){
      match(x['com1'],tax$PRIMARY_COM_NAME) -> i
      regexpr(" [(]", tax$SCI_NAME[i])[[1]][1]-1 -> p
      regexpr(" x ", tax$SCI_NAME[i])[[1]][1]-1 -> q
      #ifelse(length(p)>0,p,0)
      #ifelse(length(q)>0,q,0)
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
  
  #Add GET_MASS_CODE for com2 to the data frame (will be used to match for masses)
  sapply(d$com2,function(x){tax$SPECIES_CODE[which(tax$PRIMARY_COM_NAME==x)]},USE.NAMES=F) -> d$GET_MASS_CODE
  
  read.csv("masses_tax2016.csv",header=T,stringsAsFactors=F,na.strings=c("","NA")) -> m
  
  #Add mass & masscomm to data frame, using SPECIES_CODE
  sapply(d$GET_MASS_CODE,function(x){m[which(m$SPECIES_CODE==x),'mass']}) -> d$mass
  sapply(d$GET_MASS_CODE,function(x){m[which(m$SPECIES_CODE==x),'masscomm']}) -> d$masscomm
  "" -> d$masscomm[which(is.na(d$masscomm))]
  
  #Get subtotal per species
  d$count * d$mass -> d$subtotal
  
  #Get percent breakdown by mass
  #d$subtotal/sum(na.omit(d$subtotal))*100 -> d$pm
  d$subtotal/sum(na.omit(d$subtotal)) -> d$pm
  #sprintf("%1.2f%%", d$pm) -> d$pm
  
  #Get percent breakdown by count
  #d$count/sum(na.omit(d$count))*100 -> d$pc
  d$count/sum(na.omit(d$count)) -> d$pc
  #sprintf("%1.2f%%", d$pc) -> d$pc
  
  ## Construct final data frame
  #f <- data.frame(stringsAsFactors=F,taxon=d$com1,count=format(d$count,big.mark=","),count.p=d$pc,unit.g=format(round(d$mass,1),trim=F,big.mark=",",nsmall=1,drop0trailing=T),subtotal.g=format(round(d$subtotal,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),mass.p=d$pm,note=d$masscomm)
  f <- data.frame(stringsAsFactors=F,i=as.numeric(rownames(d)),taxon=d$com1,count=d$count,count.p=d$pc,unit.g=d$mass,subtotal.g=d$subtotal,mass.p=d$pm,note=d$masscomm)
  
  
  ## Create named list of desired objects here, and extract from result
  sum(d$count) -> ni
  sum(na.omit(d$subtotal)) -> tm
  length(d$com1) -> ns
  
  return(list('title'=title,'URL'=URL,'tm'=tm,'ni'=ni,'ns'=ns,'ifX'=ifX,'f'=f))
}

ui <- fluidPage(
  titlePanel("eBird Biomass Calculator"),
  fluidRow(
    column(3, wellPanel(
      textInput("text", "eBird Checklist URL:", ""),
      actionButton("do", "Calculate"),
      tags$br(),
      tags$br(),
      p(tags$a(href="https://github.com/slager/eBirdmass", "R code", target="_blank")," by ",tags$a(href="https://twitter.com/dlslager", "Dave Slager", target="_blank"))
    )),
    column(6,
           h5(htmlOutput("mySite")),
           DT::dataTableOutput("table"),
           h5(htmlOutput("coerceX")),
           h5(htmlOutput("nins")),
           h5(htmlOutput("metric")),
           h5(htmlOutput("english")),
           h5(htmlOutput("pennies")),
           #verbatimTextOutput("text"),
           p("Masses",tags$a(href="https://www.amazon.com/Sibley-Guide-Birds-David-Allen/dp/0679451226","(Sibley 2000)", target="_blank"),"transcribed by Sean Fitzgerald")
    )
  )
)

server <- function(input, output){
    
  eventReactive(input$do,{maketable(input$text)}) -> maketablereactive

  output$table <- DT::renderDataTable(DT::datatable({maketablereactive()[['f']]},
                                                    rownames=FALSE,
                                                    colnames=c("","Taxon","Count","%Count","Unit(g)","Subtotal(g)","%Mass","Note"),
                                                    options=list(paging=FALSE,searching=FALSE,info=FALSE)
                                                    ) %>%
                                        formatCurrency(columns=c(3,6),currency="",interval=3,mark=",",digits=0) %>%
                                        formatCurrency(columns=5,currency="",interval=3,mark=",",digits=1) %>%
                                        formatPercentage(c(4,7), 2)
                                      )
      
  #output$table <- renderTable({maketablereactive()[['f']]},align='lrrrrrl',include.rownames=FALSE)
  
  output$mySite <- renderUI({
    tags$a(href = maketablereactive()[['URL']], maketablereactive()[['title']], target="_blank")})

  output$coerceX <- renderUI({
    tags$p(style="color:red",ifelse(maketablereactive()[['ifX']],"Xs were coerced to 1s",""))
  })
  
  output$nins <- renderUI({
    tags$p(paste0(format(maketablereactive()[['ni']],big.mark=",")," individuals"),
           tags$br(),
           paste0(maketablereactive()[['ns']]," taxa"))
  })
  
  output$metric <- renderUI({
    tags$p(paste0(format(round(maketablereactive()[['tm']],1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," grams"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]/1000,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," kg"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]/1e6,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," metric tonnes"))
  })
  
  output$english <- renderUI({
    tags$p(paste0(format(round(maketablereactive()[['tm']]*0.035274,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," oz"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]*0.00220462,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," lb"),
           tags$br(),
           paste0(format(round(maketablereactive()[['tm']]*1.10231e-6,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T)," tons"))
  })
  
  output$pennies <- renderUI({
    tags$p(paste0("$",format(round(maketablereactive()[['tm']]/2.5*.01, 2), nsmall = 2, big.mark=",")," in pennies"))
  })
  
  output$text <- renderText({
    paste0(
          ifelse(maketablereactive()[['ifX']],"Xs were coerced to 1s\n\n",""),
          format(maketablereactive()[['ni']],big.mark=","),
          " individuals\n",
          maketablereactive()[['ns']],
          " taxa\n\n",
          format(round(maketablereactive()[['tm']],1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),
          " grams\n",
          format(round(maketablereactive()[['tm']]/1000,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),
          " kg\n",
          format(round(maketablereactive()[['tm']]/1e6,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),
          " metric tonnes\n\n",
          format(round(maketablereactive()[['tm']]*0.035274,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),
          " oz\n",    
          format(round(maketablereactive()[['tm']]*0.00220462,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),
          " lb\n",    
          format(round(maketablereactive()[['tm']]*1.10231e-6,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),
          " tons\n\n",
          "$",format(round(maketablereactive()[['tm']]/2.5*.01, 2), nsmall = 2, big.mark=","),
          " in pennies"
          )
  })
  
}


shinyApp(ui, server)