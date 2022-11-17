mk_stanmodel <- function(lvmodel = c("2pl","grm","gpcm")) {

  lvmodel　<- tolower(lvmodel)

  if(lvmodel == "2pl") {
    stanmodel <- ""

  } else if(lvmodel == "grm") {
    stanmodel <- ""

  } else if(lvmodel == "gpcm") {
    stanmodel <- ""

  } else {
    stop("not supported")
  }
  stanmodel
}
