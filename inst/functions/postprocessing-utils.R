check_all_diagnostics <- function(fit, outdir) {
  check_n_eff(fit)
  check_rhat(fit)
  n_div <- check_div(fit)
  n_treedepth <- check_treedepth(fit,15)
  check_energy(fit)
  
  cat('\nn_div',n_div, '\n' )
  cat('\nn_treedepth',n_treedepth, '\n' )
  saveRDS(list(n_div,n_treedepth), file=paste0(outdir,'-diagnostics.rds'))
}

check_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                n, N, max_depth, 100 * n / N))
  if (n > 0)
    print('  Run again with max_depth set to a larger value to avoid saturation')
  return(n)
}

check_energy <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)*2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('E-BFMI indicated no pathological behavior')
  else
    print('  E-BFMI below 0.2 indicates you may need to reparameterize your model')
}

check_n_eff <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  # remove na neff
  fit_summary <- fit_summary[!is.na(fit_summary[,5]),]
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(rstan::extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      print(sprintf('n_eff / iter for parameter %s is %s!',
                    rownames(fit_summary)[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('n_eff / iter looks reasonable for all parameters')
  else
    print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
}

check_div <- function(fit) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                n, N, 100 * n / N))
  if (n > 0)
    print('  Try running with larger adapt_delta to remove the divergences')
  # return iterations ended with a divergence
  return(n)
}

check_rhat <- function(fit) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  # remove na neff
  fit_summary <- fit_summary[!is.na(fit_summary[,6]),]
  N <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      print(sprintf('Rhat for parameter %s is %s!',
                    rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('Rhat looks reasonable for all parameters')
  else
    print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
}


#' Cartesian coordinates per facet-panel
#'
#' This function mimics the behavior of [ggplot2::coord_cartesian()],
#' while supporting per-panel limits when faceted.
#'
#' @details
#'
#' A 'panel_limits' data frame may contain:
#'
#' - zero or more faceting variables, all of which must be found
#'   within the grob's 'layout' (i.e., defined by
#'   [ggplot2::facet_grid()] or [ggplot2::facet_wrap()];
#' 
#' - zero or more of 'xmin', 'xmax', 'ymin', and 'ymax', where missing
#'   columns and 'NA' values within columns will default to ggplot2's
#'   normal min/max determination;
#'
#' - each panel in the plot must match no more than one row in
#'   'panel_limits';
#'
#' - each row may match more than one panel, such as when some
#'   faceting variables are not included (in 'panel_limits');
#'
#' - if no faceting variables are included, then 'panel_limits' must
#'   be at most one row (in which case it effectively falls back to
#'   [ggplot2::coord_cartesian()] behavior).
#'
#' It is an error if:
#'
#' - a panel is matched by more than one row (no matches is okay);
#' 
#' - a faceting variable in 'panel_limits' is not found within the
#'   faceted layout.
#'
#' @section Thanks:
#'
#' - burchill (github) and the original version;
#'   https://gist.github.com/burchill/d780d3e8663ad15bcbda7869394a348a
#' 
#' - Z.Lin (stackoverflow) for helping me through some of the
#'   initial errors; https://stackoverflow.com/a/63556918
#'
#' - teunbrand (github and stackoverflow), possible future extension
#'   of the non-list-index version; https://github.com/teunbrand/ggh4x
#'
#' @examples
#' \dontrun{
#' 
#' library(dplyr)
#' library(tidyr)
#' library(ggplot2)
#'
#' testdata <- tibble(
#'   x = rep(1:100, 2),
#'   y = rep(sin(seq(0,2*pi,length.out=100)), 2)
#' ) %>%
#'   mutate(y1 = y - 0.3, y2 = y + 0.3) %>%
#'   tidyr::crossing(
#'     tidyr::expand_grid(facet1 = c("aa", "bb"), facet2 = c("11", "22"))
#'   )
#' 
#' gg <- ggplot(testdata, aes(x, y)) +
#'   geom_ribbon(aes(ymin = y1, ymax = y2), fill = "#ff8888aa") +
#'   geom_path(color = "red", size = 1) +
#'   facet_wrap(facet1 + facet2 ~ ., scales = "free")
#' gg
#'
#' # single-panel change,
#' gg + coord_cartesian_panels(
#'   panel_limits = tribble(
#'     ~facet1, ~facet2, ~ymin, ~ymax
#'   , "aa"   , "22"   , -0.75,   0.5
#'   )
#' )
#'
#' # subset of facet variables, optionally tribble-style
#' gg + coord_cartesian_panels(
#'     ~facet2, ~ymin, ~ymax
#'   , "22"   , -0.75,   0.5
#' )
#' 
#' # use of 'NA' for default limits
#' gg + coord_cartesian_panels(
#'   , "aa"   ,    "11", -0.75,   0.5
#'   , "bb"   ,    "22",    NA,   0.5
#' )
#' 
#' }
#' 
#' @param panel_limits 'data.frame' with faceting variables and
#'   limiting variables, see 'Details'
#' @param expand,default,clip as defined/used in
#'   [ggplot2::coord_cartesian()]
#' @export
#' @md
coord_cartesian_panels <- function(..., panel_limits = NULL,
                                   expand = TRUE, default = FALSE, clip = "on") {
  if (is.null(panel_limits)) panel_limits <- tibble::tribble(...)
  ggplot2::ggproto(NULL, UniquePanelCoords,
                   panel_limits = panel_limits,
                   expand = expand, default = default, clip = clip)
}

UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,
  
  num_of_panels = 1,
  panel_counter = 1,
  layout = NULL,
  
  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    self$layout <- layout # store for later
    layout
  },
  
  setup_panel_params =  function(self, scale_x, scale_y, params = list()) {
    train_cartesian <- function(scale, limits, name, given_range = c(NA, NA)) {
      if (anyNA(given_range)) {
        expansion <- ggplot2:::default_expansion(scale, expand = self$expand)
        range <- ggplot2:::expand_limits_scale(scale, expansion, coord_limits = limits)
        isna <- is.na(given_range)
        given_range[isna] <- range[isna]
      }
      out <- list(
        ggplot2:::view_scale_primary(scale, limits, given_range),
        sec = ggplot2:::view_scale_secondary(scale, limits, given_range),
        arrange = scale$axis_order(),
        range = given_range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }
    
    this_layout <- self$layout[ self$panel_counter,, drop = FALSE ]
    self$panel_counter <- 
      if (self$panel_counter < self$num_of_panels) {
        self$panel_counter + 1
      } else 1
    
    # determine merge column names by removing all "standard" names
    layout_names <- setdiff(names(this_layout),
                            c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y"))
    limits_names <- setdiff(names(self$panel_limits),
                            c("xmin", "xmax", "ymin", "ymax"))
    
    limits_extras <- setdiff(limits_names, layout_names)
    if (length(limits_extras) > 0) {
      stop("facet names in 'panel_limits' not found in 'layout': ",
           paste(sQuote(limits_extras), collapse = ","))
    } else if (length(limits_names) == 0 && NROW(self$panel_limits) == 1) {
      # no panels in 'panel_limits'
      this_panel_limits <- cbind(this_layout, self$panel_limits)
    } else {
      this_panel_limits <- merge(this_layout, self$panel_limits, all.x = TRUE, by = limits_names)
    }
    
    if (isTRUE(NROW(this_panel_limits) > 1)) {
      stop("multiple matches for current panel in 'panel_limits'")
    }
    
    # add missing min/max columns, default to "no override" (NA)
    this_panel_limits[, setdiff(c("xmin", "xmax", "ymin", "ymax"),
                                names(this_panel_limits)) ] <- NA
    
    c(train_cartesian(scale_x, self$limits$x, "x",
                      unlist(this_panel_limits[, c("xmin", "xmax"), drop = TRUE])),
      train_cartesian(scale_y, self$limits$y, "y",
                      unlist(this_panel_limits[, c("ymin", "ymax"), drop = TRUE])))
  }
)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
