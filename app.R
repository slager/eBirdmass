library(shiny)
library(rvest)

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
  checklist <- paste(capture.output(html(URL), file=NULL))
  
  #Grep and trim the species
  grep('class="se-name"',checklist) -> i #grep lines of the common names
  checklist[i] -> species
  gsub("\t","",species) -> species #trim tabs
  sub('<h5 class=\"se-name\">',"",species) -> species #trim left tags
  sub('</h5>',"",species) -> species #trim right
  #Grep and trim the counts
  grep('class="se-count"',checklist) -> j #grep lines of the counts
  checklist[j] -> count
  sub('<th><h5 class=\"se-count\">',"",count) -> count #trim left
  sub('</h5></th>',"",count) -> count #trim right
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
  d$subtotal/sum(na.omit(d$subtotal))*100 -> d$pm
  sprintf("%1.2f%%", d$pm) -> d$pm
  
  #Get percent breakdown by count
  d$count/sum(na.omit(d$count))*100 -> d$pc
  sprintf("%1.2f%%", d$pc) -> d$pc
    
  ## Construct final data frame
  f <- data.frame(stringsAsFactors=F,taxon=d$com1,count=format(d$count,big.mark=","),count.p=d$pc,unit.g=format(round(d$mass,1),trim=F,big.mark=",",nsmall=1,drop0trailing=T),subtotal.g=format(round(d$subtotal,1),trim=T,big.mark=",",nsmall=1,drop0trailing=T),mass.p=d$pm,note=d$masscomm)
  
  ## Create named list of desired objects here, and extract from result
  sum(d$count) -> ni
  sum(na.omit(d$subtotal)) -> tm
  length(d$com1) -> ns
  
  return(list('tm'=tm,'ni'=ni,'ns'=ns,'ifX'=ifX,'f'=f))
}

ui <- fluidPage(
  titlePanel("eBird biomass calculator"),
  fluidRow(
    column(3, wellPanel(
      textInput("text", "eBird checklist URL:", ""),
      actionButton("do", "Calculate"),
      tags$br(),
      tags$br(),
      p(tags$a(href="https://github.com/slager/eBirdmass", "R code")," by ",tags$a(href="https://twitter.com/dlslager", "Dave Slager"))
    )),
    column(6,
           htmlOutput("mySite"),
           tableOutput('table'),
           verbatimTextOutput("text"),
           p("Masses",tags$a(href="https://www.amazon.com/Sibley-Guide-Birds-David-Allen/dp/0679451226","(Sibley 2000)"),"transcribed by Sean Fitzgerald")
    )
  )
)

server <- function(input, output){
    
  eventReactive(input$do,{maketable(input$text)}) -> maketablereactive
    
  output$table <- renderTable({maketablereactive()[['f']]},align='rlrrrrrl',include.rownames=FALSE)
  
  output$mySite <- renderUI({
    tags$a(href = input$text, "View checklist on eBird.org")})
      
  output$text <- renderText({
    paste0(
          ifelse(maketablereactive()[['ifX']],"Xs were coerced to 1s!\n\n",""),
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