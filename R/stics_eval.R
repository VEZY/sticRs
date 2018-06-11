#' Evaluate STICS version
#'
#' @description Wrapper function to run STICS simulations and evaluate
#'  outputs using observations and/or other model versions.
#'
#' @param dir.orig   Path to the directory from which to copy the simulation files. If
#'                   \code{NULL} (the default), uses the package dummy USM.
#' @param dir.targ   Path to the target directory for evaluation. Created if missing.
#' @param stics      STICS executable path named list.
#' @param obs_name   A vector of observation file name(s). It must have the form
#'                   \code{c(Dominant,Dominated)} for mixed crops.
#'                   See \code{\link{read_obs}} \code{filename} parameter for more details.
#' @param Out_var    The variables needed as output
#' @param plot_it    Boolean. Do the plot as to be pinted ?
#' @param Parallel   Are the simulations to be executed in parallel ?
#' @param mixed      (optional) Is the simulation made on mixed species (boolean)
#' @param Title      A title for the evaluation. This will be used on the ggplot object.
#'
#' @details The names in the \code{stics} parameter list are used for reference in the outputs
#' (data.frames and plots). If no names are provided, the function give a dummy name for each.
#' The function run STICS executable in parrallel using all cores from the
#' machine but one.
#'
#' @return A list of two. The first object is a list of the outputs of
#'  \code{\link{eval_output}}, which is called for each object in \code{stics}.
#'  The second is a summary of the simulations outputs as a ggplot object.
#'
#' @seealso \code{\link{eval_output}}, \code{\link{import_usm}},
#'  \code{\link{run_stics}}, and \code{\link{set_out_var}}
#'
#' @examples
#'\dontrun{
#' library(sticRs)
#' stics_eval()
#'}
#'
#' @export
stics_eval= function(dir.orig=NULL, dir.targ= getwd(),stics, obs_name= NULL,
                     Out_var=NULL, plot_it=T,Parallel=T,mixed= NULL,Title=NULL){
  if(is.list(stics)){
    usm_name= names(stics)
  }else{
    stop("stics parameter must be a list")
  }

  if(Parallel){
    NbCores= parallel::detectCores()-1
    cl= parallel::makeCluster(min(NbCores,length(stics)))
    parallel::clusterExport(cl=cl,
                  varlist=c("dir.orig","dir.targ","usm_name","stics",
                            "obs_name","Out_var","import_usm","set_out_var",
                            "run_stics","eval_output"),
                  envir=environment())
    outputs=
      parallel::parLapply(cl,seq_along(stics),
                function(x,dir.orig,dir.targ,usm_name,stics,obs_name){
                  USM_path= file.path(dir.targ,usm_name[x])
                  import_usm(dir.orig = dir.orig, dir.targ = dir.targ,
                             usm_name = usm_name[x], overwrite = T,
                             stics = stics[[x]])
                  set_out_var(filepath= file.path(USM_path,"var.mod"),
                              vars=Out_var, app=F)
                  run_stics(dirpath = USM_path)
                  output= eval_output(dirpath= USM_path, obs_name= obs_name)
                  output
                },dir.orig,dir.targ,usm_name,stics,obs_name)
    parallel::stopCluster(cl)
  }else{
    outputs=
      lapply(seq_along(stics),
             function(x){
               USM_path= file.path(dir.targ,usm_name[x])
               import_usm(dir.orig = dir.orig, dir.targ = dir.targ,
                          usm_name = usm_name[x], overwrite = T,
                          stics = stics[[x]])
               set_out_var(filepath= file.path(USM_path,"var.mod"),
                           vars=Out_var, app=F)
               run_stics(dirpath = USM_path)
               output= eval_output(dirpath= USM_path, obs_name= obs_name)
               output
             })
  }
  names(outputs)= usm_name
  outputs[["plot_it"]]= plot_it
  outputs[["Vars"]]= Out_var
  outputs[["Title"]]= Title

  gg_stics= do.call(plot_output,outputs)
  return(list(outputs= outputs[-grep("plot_it|Vars|Title",names(outputs))],
              gg_object= gg_stics))
}