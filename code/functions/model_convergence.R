# =============================================================================
# code/functions/model_convergence.R
#
# Helper functions to check convergence of fitted GAM/GLM models.
# Inspects model convergence flags in mgcv::gam fits to prevent silent fitting
# failures.
# =============================================================================

#' Check mgcv::gam Model Convergence
#'
#' Evaluates whether a fitted gam model converged successfully, inspecting both
#' the inner algorithm status and outer iteration info.
#'
#' @param model_obj A fitted gam model from mgcv.
#' @param model_name Character. Optional name of the model for logging.
#' @param raise_warning Logical. If TRUE, raises a warning if model did not converge.
#'
#' @return A list containing:
#'   \describe{
#'     \item{converged}{Logical. TRUE if converged, FALSE otherwise.}
#'     \item{message}{Character. Descriptive convergence status message.}
#'   }
#' @export
check_model_convergence <- function(model_obj, model_name = "Model", raise_warning = TRUE) {
  if (is.null(model_obj)) {
    return(list(converged = FALSE, message = "Model object is NULL"))
  }
  
  converged <- TRUE
  msg <- "Model converged successfully."
  
  # Check inner convergence
  if (!is.null(model_obj$converged) && !model_obj$converged) {
    converged <- FALSE
    msg <- "Inner algorithm did not converge."
  }
  
  # Check outer convergence info if present
  if (converged && !is.null(model_obj$outer.info)) {
    outer_conv <- model_obj$outer.info$conv
    if (!is.null(outer_conv) && is.character(outer_conv)) {
      if (grepl("non|fail|warn|limit", tolower(outer_conv)) ||
          !grepl("conv|full|ok", tolower(outer_conv))) {
        converged <- FALSE
        msg <- paste("Outer optimization warning:", outer_conv)
      }
    }
  }
  
  if (!converged && raise_warning) {
    warning(sprintf("Convergence failure in model '%s': %s", model_name, msg))
  }
  
  return(list(converged = converged, message = msg))
}
