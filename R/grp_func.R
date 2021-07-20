#' Group ID Function - Double Sample Approach
#' Implementation described in est_fncn()
#'
#' @param LT Tree List
#' @param id ID
#' @export
#' @keywords Tree Subgroup
#'

grp_func <- function(LT, id){
  if(LT$nodetype == "<leaf>"){
    grp.ids <- LT$index.TEST
  }else{
    k=1
    str <- "LT$"
    while(k < 100){
      if(tryCatch(id %in% eval(parse(text=paste(str,"left$index.TEST",sep=""))),
                  error = function(i){return(FALSE)})){
        grp.ids <- eval(parse(text=paste(str,"left$index.TEST",sep="")))
        str <- paste(str,"left$",sep="")
        k=k+1
      }else if(tryCatch(id %in% eval(parse(text=paste(str,"right$index.TEST",sep=""))),
                        error=function(i){return(FALSE)})){
        grp.ids <- eval(parse(text=paste(str,"right$index.TEST",sep="")))
        str <- paste(str,"right$",sep="")
        k=k+1
      }else{
        break
      }
    }
  }
  return(grp.ids)
}
