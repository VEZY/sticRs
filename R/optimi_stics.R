#' STICS parameter optimization
#'
#' @description Optimize STICS parameter values according to measurements
#'
#' @param dir.orig   Path to the directory from which to copy the simulation files. If
#'                   \code{NULL} (the default), uses the package dummy USM.
#' @param dir.targ   Path to the target directory for evaluation. Created if missing.
#' @param stics      STICS executable path
#' @param obs_name   A vector of observation file name(s). It must have the form
#'                   \code{c(Dominant,Dominated)} for mixed crops.
#'                   See \code{\link{read_obs}} \code{filename} parameter for more details.
#' @param Parameters A list of list of starting, min and max values for each parameters, named after
#'                   them (see details and example)
#' @param Vars       Output variables on which the optimization is performed
#' @param weight   The weight used for each variable (see details)
#' @param method     The optimization method to use, see \pkg{dfoptim} package. For the moment, only \code{\link[dfoptim]{nmkb}}
#' @param Plant      The plant (\emph{i.e.} Principal or associated) for which the parameters
#'                   will be set (only for plant or technical parameters in mixed crop simulations)
#'                   Set to \code{NULL} if using STICS in sole crop
#' @param ...        Further parameters passed to the optimization function called
#'                   (see \pkg{dfoptim} package)
#'
#' @details The function uses the \pkg{dfoptim} package functions under the hood. Currently only the Nelder-Mead algorithm
#'  is implemented.
#'  The `Parameters` argument should take the form of a list of arguments for each parameter, named after the parameter
#'  of interest (see example). If the start is `NULL`, then the mean value between the min and max values is taken.
#'  If weight is not provided by the user, the selection criteria is computed using the equation
#' 5 from Wallach et al. (2011). If they are provided, the equation 6 is used instead.
#'
#' @references Wallach, D., Buis, S., Lecharpentier, P., Bourges, J., Clastre, P., Launay, M., … Justes, E. (2011).
#'  A package of parameter estimation methods and implementation for the STICS crop-soil model. Environmental Modelling & Software, 26(4), 386–394. doi:10.1016/j.envsoft.2010.09.004
#'
#'
#' @return A list of three :
#' \itemize{
#'   \item gg_objects: A list of ggplot objects to plot the final STICS simulation with optimized parameter values
#'   compared to original parameter values.
#'   \item values: A list of the optimized parameter values.
#' }
#'
#' @importFrom dfoptim nmkb
#' @importFrom dplyr group_by summarise summarise_all select
#' @importFrom magrittr "%<>%"
#'
#' @seealso \code{\link[dfoptim]{nmkb}}
#'
#' @examples
#'\dontrun{
#' library(sticRs)
#'
#' Parameters= data.frame(parameter= c('hautK1','hautK2'),
#'                        start= c(0.2,0.2),
#'                        min= c(0,0),
#'                        max= c(1,1))
#'
#' optimi_stics(dir.orig = "0-DATA/dummy/Year_2005_2006/IC_Wheat_Pea",
#'              dir.targ = "2-Simulations/param_optim",
#'              stics = "0-DATA/stics_executable/19-new/Stics.exe",
#'              obs_name = c("6_IC_Wheat_N0.obs","6_IC_Pea_N0.obs"),
#'              Parameters= Parameters,
#'              Vars = c('lai(n)','masec(n)','hauteur'),
#'              method= "nmkb",Plant=1)
#'}
#'
#' @export
#'
optimi_stics= function(dir.orig, dir.targ=getwd(),stics,obs_name,Parameters,
                       Vars,weight=NULL,method=c("nmkb"),Plant=1,...){
  .= NULL # to avoid CRAN checks errors

  method= match.arg(method,c("nmkb")) # add new methods here

  Vars_R= gsub("\\(","_",Vars)%>%gsub("\\)","",.)

  if(!is.null(weight)){
    if(length(weight)!=length(Vars_R)){
      stop("weight length should be the same as the Vars length")
    }
    weight= data.frame(variable= Vars_R, weight= weight, stringsAsFactors = FALSE)
  }

  if(!is.data.frame(Parameters)){
    stop("Parameters should be a data.frame")
  }else if(!all(colnames(Parameters)%in%c("parameter","start","min","max"))){
    stop("Parameters should be a data.frame with columns name: parameter, start, min, max")
  }else{
    null_start= lapply(Parameters$start,is.null)%>%unlist
    if(any(null_start)){
      Parameters$start[null_start]= (Parameters$max[null_start]+Parameters$min[null_start])/2
    }
  }

  # NB: we don't use stics_eval dircectly because we want to copy the usm only once for
  # performance.

  usm_name= paste0("optim")
  USM_path= file.path(dir.targ,usm_name)
  import_usm(dir.orig = dir.orig, dir.targ = dir.targ,
             usm_name = usm_name, overwrite = T,
             stics = stics)
  set_out_var(filepath= file.path(USM_path,"var.mod"),
              vars=Vars, add=F)
  run_stics(dirpath = USM_path)
  output_orig= eval_output(dirpath= USM_path, obs_name= obs_name)

  opti= dfoptim::nmkb(fn= stics_eval_opti,
                      par= Parameters$start,
                      lower= Parameters$min,
                      upper= Parameters$max,
                      USM_path= USM_path,
                      param= as.character(Parameters$parameter),
                      Plant= Plant,
                      weight= weight,
                      obs_name= obs_name)

  # Import the last simulation output:
  output_opti= eval_output(dirpath= USM_path, obs_name= obs_name)
  gg_output= plot_output(original= output_orig,
                         optimized= output_opti, plot_it = FALSE)

  unlink(x = dir.targ, recursive = T, force = T)

  invisible(list(gg_objects= gg_output, opti_output= opti,last_sim= output_opti))
}




#' Objective function to optimize
#'
#' @description Function evaluated for parameter optimization. This function is only provided
#' for informative value, but should not be used by the user. It is called by \code{\link{optimi_stics}}.
#'
#' @param x        The starting parameter values
#' @param USM_path The path to the USM
#' @param obs_name The observation file name
#' @param param    The parameter names (in same order than x)
#' @param weight   The weight used for each variable (see details)
#' @param Plant    The plant evaluated (for intercrop)
#'
#' @details If weight is not provided by the user, the selection criteria is computed using the equation
#' 5 from Wallach et al. (2011). If they are provided, the equation 6 is used instead.
#'
#' @references Wallach, D., Buis, S., Lecharpentier, P., Bourges, J., Clastre, P., Launay, M., … Justes, E. (2011).
#'  A package of parameter estimation methods and implementation for the STICS crop-soil model. Environmental Modelling & Software, 26(4), 386–394. doi:10.1016/j.envsoft.2010.09.004
#'
#' @importFrom rlang .data
#'
#' @return The weighted product of squares (selection criteria)
#' @export
#'
stics_eval_opti= function(x,USM_path,obs_name,param,weight=NULL,Plant=1){
  names(x)= param
  lapply(param, function(pa){
    set_param(dirpath = USM_path, param = pa,
              value = x[pa], plant = Plant)
  })
  run_stics(dirpath = USM_path)
  out_stats=
    stati_stics(USM_path, obs_name= obs_name)%>%
    dplyr::ungroup()

  # Selecting only the plant the user need:
  if(length(unique(out_stats$Dominance))>1){
    out_stats= out_stats[out_stats$Dominance==ifelse(Plant==1,"Principal","Associated"),]
  }

  if(is.null(weight)){
    crit=
      out_stats%>%
      dplyr::mutate(crit= (.data$SS_res/.data$n_obs)^(.data$n_obs/2))%>%
      dplyr::summarise(crit= prod(.data$crit))
  }else{
    crit=
      dplyr::left_join(weight,out_stats%>%
                         dplyr::mutate(variable= as.character(.data$variable)),
                       by= "variable")%>%
      dplyr::mutate(crit= .data$SS_res*.data$weight)%>%
      dplyr::summarise(crit= sum(.data$crit))
  }
  crit$crit
}

