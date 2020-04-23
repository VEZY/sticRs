#' Run STICS
#'
#' @description Helper function to start a simulation using
#' STICS executable
#'
#' @param dirpath USM directory path
#' @param exe_name (optional) The executable name.
#'
#' @return The function prints the STICS output to the console and returns
#'         `TRUE` if STICS ran successfully, or an error if any problem occured
#'
#' @examples
#'\dontrun{
#' library(sticRs)
#' run_stics()
#'}
#' @export
#'
run_stics=function(dirpath=getwd(),exe_name="stics"){
  wd= getwd()
  on.exit(setwd(wd))
  setwd(dirpath)
  out= system2(exe_name)
  if(out==0){
    TRUE
  }else{
    stop("STICS call ",crayon::red("failed"))
  }
}
