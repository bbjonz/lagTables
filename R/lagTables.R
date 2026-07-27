#' @importFrom stats anova as.formula chisq.test deviance df.residual drop1 ftable glm na.omit pchisq poisson xtabs
#' @importFrom utils data
NULL

#these are columns/list-elements created via dplyr/tidyr non-standard
#evaluation inside the functions below (not undefined globals) -- listed
#here so R CMD check doesn't flag them
utils::globalVariables(c(
  ".", ".data", ":=", "group", "tempcol", "Previous Unit(s)", "Var1", "Var2", "Freq",
  "freq", "data", "tidied", "time12", "model.fit", "model.stats", "fit.glanced",
  "deviance", "df.residual", "model.p", "devp.mod", "Term"
))

#' Add lag columns of a variable to a data frame
#'
#' Internal helper used by \code{trprobs}/\code{lagmodels}/\code{shmodels}
#' to build the lag1..lagN columns needed for lag sequential analysis. If
#' \code{d} is grouped (via \code{dplyr::group_by}), lags correctly reset
#' at each group boundary instead of bleeding across groups -- this
#' matters whenever \code{laggroup} identifies separate sequences (e.g.
#' separate meetings) that should not be lagged into one another.
#'
#' @param d a (possibly grouped) data frame
#' @param col name of the column to lag (character)
#' @param lagnum number of lags to create
#' @return \code{d} with lag1..lagN columns added (zero-padded names,
#'   e.g. \code{lag01}, \code{lag02} for \code{lagnum >= 10})
#' @noRd
add_lags <- function(d, col, lagnum) {
  lag_names <- paste0("lag", formatC(seq_len(lagnum), width = nchar(lagnum), flag = "0"))
  for (i in seq_len(lagnum)) {
    d <- dplyr::mutate(d, "{lag_names[i]}" := dplyr::lag(.data[[col]], i))
  }
  d
}

#' Format a numeric data frame to a fixed number of decimals
#'
#' Base R replacement for \code{DescTools::Format(x, digits = digits)}.
#' Appropriate for statistics that can exceed 1 in magnitude (expected
#' counts, standardized residuals, deviances) -- these should keep a
#' leading zero when the value is less than 1 (e.g. "0.23"). See
#' \code{fmt_bounded} for statistics naturally bounded to [-1, 1].
#'
#' @param x a data frame of numeric columns
#' @param digits number of decimal places
#' @return a data frame of character strings, same dimensions as \code{x}
#' @noRd
fmt_fixed <- function(x, digits = 2) {
  out <- lapply(x, function(col) sprintf(paste0("%.", digits, "f"), col))
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  names(out) <- names(x)
  rownames(out) <- rownames(x)
  out
}

#' Format a numeric data frame to a fixed number of decimals, dropping
#' the leading zero
#'
#' For statistics naturally bounded to [-1, 1] -- correlations,
#' proportions, transition probabilities -- where ".23" is the usual
#' convention rather than "0.23". Not appropriate for statistics that
#' can exceed 1 in magnitude; use \code{fmt_fixed} for those.
#'
#' @param x a data frame of numeric columns
#' @param digits number of decimal places
#' @return a data frame of character strings, same dimensions as \code{x}
#' @noRd
fmt_bounded <- function(x, digits = 2) {
  out <- fmt_fixed(x, digits)
  out[] <- lapply(out, function(col) sub("^(-?)0\\.", "\\1.", col))
  out
}


#' Create Transitional Probability Tables and Plots
#'
#' This function works with a vector of sequential data and creates
#' html tables for the observed and expected frequencies, as well
#'
#' If requested, it will create bar plots of standardized residuals
#' for transitional probs (Lehmann-Willenbrock, 2013)
#'
#' Finally, it will create tables and plots by group if requested
#'
#' Default is one lag but higher-order lags can be requested
#'
#' @param d dataframe that contains the sequential vector
#' @param lagvar the name of the vector that contains sequential data
#' @param laggroup vector that contains the grouping variables (default is NULL)
#' @param lagnum the number of lags to report.  Default is 1
#' @param plots request bar plots of standardized residuals. Default is 0. Set to 1 for plots.
#' @param title specify title for tables and plots.
#' @param dname provide a name for the data set for table caption.  Default is blank.
#' @return Table with tables for the observed and expected frequencies, as well as the transitional probabilities and standardized residuals.  Bar plots if requested
#' @export

trprobs <- function(d, lagvar, laggroup=NULL, lagnum=1, plots=0, dname="",
                    title="Lag Sequential Descriptive Stats") {

  options(scipen = 999, warn = -1)

  d <- as.data.frame(d, stringsAsFactors=FALSE)

  if(length(d[lagvar])==0) {stop("The variable does not exist in your data frame")}

  `%>%` <- magrittr::`%>%`

  ####IF GROUPING VARIABLE ####
  if(length(d[laggroup])!=0) {

    #lags are computed within each group so they don't bleed across
    #group boundaries (e.g. from the end of one meeting into the start
    #of the next)
    lagdat <- d %>%
      dplyr::rename(lag0 = dplyr::all_of(lagvar), group = dplyr::all_of(laggroup)) %>%
      dplyr::group_by(group) %>%
      add_lags("lag0", lagnum) %>%
      dplyr::ungroup()

    lagdat %>%
      dplyr::select(group, dplyr::contains("lag")) %>%
      ftable(xtabs(, data=.)) %>%
      as.matrix() %>%
      chisq.test(, simulate.p.value = TRUE) -> lag.tab

    obs <- as.data.frame(lag.tab$observed)
    #expect/stdres are unbounded (can exceed 1) -- keep the leading zero.
    #tr is a transition probability, always in [0,1] -- drop it, per the
    #usual convention for bounded/correlation-like statistics
    expect <- fmt_fixed(as.data.frame(lag.tab$expected), digits = 2)
    stdres <- fmt_fixed(as.data.frame(lag.tab$stdres), digits = 2)
    tr <- fmt_bounded(as.data.frame(obs/rowSums(obs)), digits = 2)

    mypaste <- function(x,y) paste(x, "<br>(", y, ")", sep="")

    obs.exp <- mapply(mypaste, obs, expect)

    tr.std <- mapply(mypaste, tr, stdres)

    rownames(tr.std) <- rownames(obs.exp) <- rownames(obs)

    numcodes <- nrow(obs.exp)

    print(
      data.frame(tempcol = row.names(rbind(obs.exp, tr.std)),
                 rbind(obs.exp, tr.std)) %>%
        tibble::as_tibble() %>%
        tibble::remove_rownames() %>%
        tidyr::separate(tempcol, c("Group","Previous Unit(s)"),
                        sep="_", extra="merge") %>%
        dplyr::mutate(`Previous Unit(s)` = sub("_", "&#8594;",
                                                `Previous Unit(s)`, fixed = TRUE)) %>%
        knitr::kable(caption=paste0("Lag ", lagnum, "Transition Probabilities By Group"),
                     escape = F) %>%
        kableExtra::kable_styling(bootstrap_options = c("striped", "hover"),
                                  full_width = F) %>%
        #collapse_rows() has to run *before* pack_rows() here -- calling it
        #after inserting the pack_rows() group-header rows makes kableExtra
        #recompute row groups in a way that silently blanks the label text
        #on the second (and any later) pack_rows() call
        kableExtra::collapse_rows(columns = 1, valign = "top") %>%
        kableExtra::pack_rows("Observed Frequencies\n(Expected Frequencies)", 1,
                              numcodes) %>%
        kableExtra::pack_rows("Transitional Probablities\n(Standardized Residuals)",
                              numcodes+1, numcodes*2) %>%
        kableExtra::add_header_above(c(" " = 2, "Target Unit" = ncol(obs.exp)))
    )

    #plot if requested
    if(plots>0) {
      print(
        lag.tab$stdres %>%
          as.table() %>%
          as.data.frame() %>%
          dplyr::rename(tempcol = 1, Var2 = 2) %>%
          tidyr::separate(tempcol, c("Group", "Var1"), sep="_", extra="merge") %>%
          dplyr::group_by(Var1) %>%
          ggplot2::ggplot(ggplot2::aes(x=Var2, y=Freq, fill=Var2)) +
          ggplot2::geom_bar(stat='identity', width = .4) +
          ggplot2::geom_hline(yintercept = 1.96, linetype="dashed") +
          ggplot2::geom_hline(yintercept = -1.96, linetype="dashed") +
          ggplot2::coord_flip() +
          ggplot2::facet_grid(Group ~ Var1) +
          ggplot2::scale_fill_grey() +
          ggplot2::labs(y="", x = "", fill = "Unit\nType",
                        title = paste0("Standardized Residuals For All Units"))
      )
    }

    ####IF NO GROUP VARIABLE, PROCESS DATA WITHOUT GROUPS####
  } else {

    lagdat <- d %>%
      dplyr::rename(lag0 = dplyr::all_of(lagvar)) %>%
      add_lags("lag0", lagnum)

    lagdat %>%
      dplyr::select(dplyr::contains("lag")) %>%
      ftable(xtabs(, data=.)) %>%
      as.matrix() %>%
      chisq.test(, simulate.p.value = TRUE) -> lag.tab

    obs <- as.data.frame(lag.tab$observed)
    #expect/stdres are unbounded (can exceed 1) -- keep the leading zero.
    #tr is a transition probability, always in [0,1] -- drop it, per the
    #usual convention for bounded/correlation-like statistics
    expect <- fmt_fixed(as.data.frame(lag.tab$expected), digits = 2)
    stdres <- fmt_fixed(as.data.frame(lag.tab$stdres), digits = 2)
    tr <- fmt_bounded(as.data.frame(obs/rowSums(obs)), digits = 2)

    mypaste <- function(x,y) paste(x, "<br>(", y, ")", sep="")

    obs.exp <- mapply(mypaste, obs, expect)

    tr.std <- mapply(mypaste, tr, stdres)

    rownames(tr.std) <- rownames(obs.exp) <- rownames(obs)

    numcodes <- nrow(obs.exp)

    print(
      data.frame(tempcol = row.names(rbind(obs.exp, tr.std)), rbind(obs.exp, tr.std)) %>%
        tibble::as_tibble() %>%
        tibble::remove_rownames() %>%
        dplyr::mutate(tempcol = sub("_", "&#8594;", tempcol, fixed = TRUE)) %>%
        dplyr::rename("Previous Unit(s)"=tempcol) %>%
        knitr::kable(caption=paste0("Lag ", lagnum, " Transition Probabilities"),
                     escape = F) %>%
        kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>%
        kableExtra::pack_rows("Observed Frequencies\n(Expected Frequencies)", 1, numcodes) %>%
        kableExtra::pack_rows("Transitional Probablities\n(Standardized Residuals)", numcodes+1, numcodes*2) %>%
        kableExtra::add_header_above(c(" " = 1, "Target Unit" = ncol(obs.exp)))
    )

    #plot if requested
    if(plots>0) {
      print(
    lag.tab$stdres %>%
      as.table() %>%
      as.data.frame() %>%
      dplyr::rename(Var1 = 1, Var2 = 2) %>%
      dplyr::group_by(Var1) %>%
      ggplot2::ggplot(ggplot2::aes(x=Var2, y=Freq, fill=Var2)) +
      ggplot2::geom_bar(stat='identity', width = .5) +
      ggplot2::geom_hline(yintercept = 1.96, linetype="dashed") +
      ggplot2::geom_hline(yintercept = -1.96, linetype="dashed") +
      ggplot2::coord_flip() +
      ggplot2::facet_wrap(~ Var1) +
      ggplot2::scale_fill_grey() +
      ggplot2::labs(y="Standardized Residuals", x = "", fill = "Comment\nType",
                    title = paste0("Standardized Residuals For All Comment Types"))
      )
    }

    }
  #End function
}


#' Estimates log linear models for stationarity of sequential data
#' Requires a vector of sequential data
#' Can request estimates by group
#' @param d dataframe containing sequential vector
#' @param lagcol vector in dataframe containing sequential vector
#' @param laggroup vector in dataframe containing grouping variable.  Default is empty (no groups)
#' @param lagnum number of lags. Default is 1
#' @param title caption for table
#' @return Table with log linear estimates for lag sequential data
#' @export

lagmodels <- function(d, lagcol, laggroup="", title="Log Linear Models for Stationarity", lagnum=1) {

  options(scipen = 999, warn = -1)

  if(length(lagcol)==0) {stop("The variable does not exist in your data frame")}

  #define magrittr pipe
  `%>%` <- magrittr::`%>%`

#####Group NOT Specified#####
  if(missing(laggroup)) {

    lag.counts <- d %>%
      dplyr::select(dplyr::all_of(lagcol)) %>%
      dplyr::rename(lag0 = dplyr::all_of(lagcol)) %>%
      add_lags("lag0", lagnum) %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = dplyr::n(), .groups = "drop") %>%
      as.data.frame()

    #the ncol routine removes freq from the end of the formula
    lag.form <- as.formula(paste("freq", paste(colnames(lag.counts[-(ncol(lag.counts))]),
                                               collapse=" * "), sep="~"))

    fit <- glm(lag.form, data=lag.counts, family=poisson)

    fit.out <- anova(fit, test="Chisq")

    suppressWarnings(
      broom::tidy(fit.out) %>%
        knitr::kable(digits = 2, align=c("l",rep("c",5)),
                     col.names = c("Model",
                                   "Model df",
                                   "Deviance",
                                   "Residual Deviance df",
                                   "Residual Deviance",
                                   "p value")) %>%
        kableExtra::kable_styling(full_width = F))

#####Group Specified#####
  } else {

    if(length(laggroup)==0) {stop("The group variable does not exist in your data frame")}

    #lags are computed within each laggroup level so they don't bleed
    #across group boundaries (e.g. from the end of one meeting into the
    #start of the next)
    lag.counts <- d %>%
      dplyr::select(dplyr::all_of(c(lagcol, laggroup))) %>%
      dplyr::rename(lag0 = dplyr::all_of(lagcol)) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(laggroup))) %>%
      add_lags("lag0", lagnum) %>%
      dplyr::ungroup() %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = dplyr::n(), .groups = "drop") %>%
      as.data.frame()

    lag.counts %>%
      dplyr::select(-freq, -dplyr::all_of(laggroup)) -> lag.names

    lag.form <- as.formula(paste("freq", paste(colnames(lag.names),
                                               collapse=" * "), sep="~"))
    options(knitr.kable.NA = '')
    #https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
    lag.counts %>%
      dplyr::arrange_at(dplyr::vars(dplyr::one_of(laggroup))) %>%
      na.omit() %>%
      tidyr::nest(-laggroup) %>%
      dplyr::mutate(
        tidied = purrr::map(data, ~
                              glm(lag.form, family = poisson, data = .x) %>%
                              anova(., test = "Chisq") %>%
                              broom::tidy(.))) %>%
      tidyr::unnest(tidied) %>%
      dplyr::select(-data) %>% #removes the data column

      knitr::kable(digits = 2,
                   col.names = c(laggroup,
                                 "Model",
                                 "Model df",
                                 "Deviance",
                                 "Residual Deviance df",
                                 "Residual Deviance",
                                 "p value")) |>
      kableExtra::kable_styling(full_width = F)
  }

  #end function
}


#' Tests stationarity/homogeneity across groups for sequential data
#' Requires a vector of sequential data
#' Requires a group variable
#' Tests first half of discussion against second half
#' @param d dataframe containing sequential vector
#' @param lagcol vector in dataframe containing sequential vector
#' @param laggroup vector in dataframe containing grouping variable.  Default is empty (no groups)
#' @param lagnum number of lags. Default is 1
#' @param title caption for table
#' @return Table with log linear estimates for lag sequential data
#' @export

shmodels <- function(d, lagcol, laggroup="", title="Homogeniety Tests for Log Linear Models", lagnum=1) {

  #https://data.library.virginia.edu/an-introduction-to-loglinear-models/

  options(scipen = 999, warn = -1)

  #make sure sequential data column exists
  if(length(d[[lagcol]])==0) {stop("The variable does not exist in your data frame")}

  #define magrittr pipe
  `%>%` <- magrittr::`%>%`

#####No Grouping Variable#####
  if(missing(laggroup)) {

    lag.counts <- d %>%
      dplyr::select(dplyr::all_of(lagcol)) %>%
      dplyr::mutate(time12 = ifelse(dplyr::row_number() < dplyr::n()/2, 1, 2)) %>%
      dplyr::rename(lag0 = dplyr::all_of(lagcol)) %>%
      add_lags("lag0", lagnum) %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = dplyr::n(), .groups = "drop") %>%
      as.data.frame()

    lag.counts %>%
      dplyr::select(-freq, -time12) -> lag.names

    post.form <- paste(colnames(lag.names), collapse  = ":")
    post0.form <- paste(colnames(lag.names[-length(lag.names)]), collapse  = ":")
    pre.form <- paste("freq", paste(colnames(lag.names), collapse="*"), sep="~")
    lag.form <- as.formula(paste0(pre.form, "*time12", " - ",post.form,":time12"))
    ant.form <- paste0("freq", "~", post.form, " + ", post0.form, ":time12")
    model.form <- as.formula(ant.form)

    options(knitr.kable.NA = '')
    #https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
    lag.counts %>%
      na.omit() %>%
        glm(model.form, family = poisson, data = .) -> fit.out

    lag.counts %>%
      na.omit() %>%
        glm(lag.form, family = poisson, data = .) %>%
        anova(., test = "Chisq") -> stats.out

    broom::tidy(stats.out) %>%
      knitr::kable(digits = 2, escape = F, format="html",
                   col.names = c(
                                 "Model",
                                 "Model df",
                                 "Deviance",
                                 "Residual Deviance df",
                                 "Residual Deviance",
                                 "p value")) %>%
      kableExtra::kable_styling(full_width = F) %>%
      kableExtra::add_footnote(paste0("Note: Model deviance = ",
                                      round(broom::glance(fit.out)$deviance,2),", p = ",
                                      round(pchisq(broom::glance(fit.out)$deviance,
                                             df = broom::glance(fit.out)$df.residual,
                                             lower.tail = F),2),". See Poole (2000)."))
    #ends non-group function

#####Grouping Variable#####
    } else {

      #Make sure grouping variable exists
  if(length(d[[laggroup]])==0) {stop("The group variable does not exist in your data frame")}

  #handles a weird quoting issue with multiple calls to
  #the same var in dplyr
  laggroup2 <- rlang::ensyms(laggroup)

    #time12 (first half/second half) is computed per-group, and lags are
    #also computed per-group (kept grouped through add_lags) so neither
    #the halving nor the lags bleed across group boundaries
    lag.counts <- d %>%
      dplyr::select(lagcol, laggroup) %>%
      dplyr::group_by(!!!laggroup2) %>%
      dplyr::mutate(time12 = ifelse(dplyr::row_number() < dplyr::n()/2, 1, 2)) %>%
      dplyr::rename(lag0 = dplyr::all_of(lagcol)) %>%
      add_lags("lag0", lagnum) %>%
      dplyr::ungroup() %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = dplyr::n(), .groups = "drop") %>%
      as.data.frame()

    lag.counts %>%
      dplyr::select(-freq, -dplyr::all_of(laggroup), -time12) -> lag.names

    post.form <- paste(colnames(lag.names), collapse  = ":")
    post0.form <- paste(colnames(lag.names[-length(lag.names)]), collapse  = ":")
    pre.form <- paste("freq", paste(colnames(lag.names), collapse="*"), sep="~")
    lag.form <- as.formula(paste0(pre.form, "*time12", " - ",post.form,":time12"))
    ant.form <- paste0("freq", "~", post.form, " + ", post0.form, ":time12")
    model.form <- as.formula(ant.form)

    options(knitr.kable.NA = '')
    #https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html
    lag.counts %>%
      dplyr::arrange_at(dplyr::vars(dplyr::one_of(laggroup))) %>%
      na.omit() %>%
      tidyr::nest(-laggroup) %>%
      dplyr::mutate(
        model.fit = purrr::map(data, ~ glm(model.form, family = poisson, data = .x)),
        fit.glanced = purrr::map(model.fit, broom::glance),
        model.stats = purrr::map(data, ~ glm(lag.form, family = poisson, data = .x) %>%
                          anova(., test = "Chisq") %>%
                            broom::tidy(.))) -> test.out

    #tidyr::unnest()'s old .drop=TRUE argument used to drop the other
    #list-columns automatically; it's a no-op as of tidyr 1.0 ("all
    #list-columns are now preserved"), so the leftover list-columns are
    #dropped explicitly here instead -- otherwise they survive into the
    #kable() call below and break its column-name assignment
    test.out %>%
      tidyr::unnest(model.stats) %>%
      dplyr::select(-data, -model.fit, -fit.glanced) -> stats.out

    test.out %>%
     tidyr::unnest(fit.glanced) %>%
      dplyr::mutate(model.p = pchisq(deviance, df = df.residual, lower.tail = F)) %>%
      dplyr::mutate(devp.mod = paste0(round(deviance,2),"<br>(p. = ",round(model.p),2,")")) %>%
      dplyr::select(laggroup, devp.mod) -> fit.out

    dplyr::left_join(fit.out, stats.out) %>%
      knitr::kable(digits = 2, escape = F, format="html",
                   col.names = c(laggroup,
                                 "Model Deviance",
                                 "Model",
                                 "Model df",
                                 "Deviance",
                                 "Residual Deviance df",
                                 "Residual Deviance",
                                 "p value")) %>%
      kableExtra::kable_styling(full_width = F) %>%
      kableExtra::collapse_rows(columns = 1:2, valign = "top") %>%
      kableExtra::add_footnote("Note: Model deviance is the test for stationarity as described in Poole (2000)")

    #ends group function
    }

  #end function
}


#' Tests homogeneity of sequential structure across groups
#'
#' Requires a vector of sequential data and a grouping variable (e.g.
#' separate meetings/sessions). Unlike \code{shmodels()}, which splits each
#' sequence in half to test whether structure holds constant over time,
#' \code{hommodels()} fits a single joint log-linear model across all
#' groups and tests whether the lag structure itself differs between
#' groups (the highest-order lag-by-group interaction) -- e.g. whether
#' the same transition pattern holds across several meetings that share
#' the same identified lag order. Per Poole (2000), \code{lagnum} should
#' include one lag beyond the order already identified for the sequences
#' being tested (e.g. \code{lagnum=2} for sequences with an identified
#' lag1 structure, \code{lagnum=3} for lag2), so the partial/control term
#' is present in the model. Since there is no such thing as homogeneity
#' with no groups to compare, \code{laggroup} is required (unlike the
#' other functions in this package, where it is optional).
#'
#' @param d dataframe containing sequential vector
#' @param lagcol vector in dataframe containing sequential vector
#' @param laggroup vector in dataframe containing the grouping variable
#'   whose levels are tested for homogeneity (e.g. \code{"Meeting"}).
#'   Required -- lags are also reset at each level's start, so this same
#'   variable should identify separate sequences.
#' @param lagnum number of lags. Default is 1
#' @param title caption for table
#' @return Table with the log-linear homogeneity test (sequential
#'   \code{anova}) for the lag structure across \code{laggroup} levels
#' @export

hommodels <- function(d, lagcol, laggroup, lagnum=1, title="Homogeneity Tests for Log Linear Models") {

  options(scipen = 999, warn = -1)

  #make sure sequential data column exists
  if(length(d[[lagcol]])==0) {stop("The variable does not exist in your data frame")}

  #homogeneity is a comparison across groups, so a group variable is
  #always required (unlike trprobs/lagmodels/shmodels, where it's optional)
  if(missing(laggroup) || length(d[[laggroup]])==0) {
    stop("A grouping variable (laggroup) is required to test homogeneity across groups")
  }

  #define magrittr pipe
  `%>%` <- magrittr::`%>%`

  #handles a weird quoting issue with multiple calls to
  #the same var in dplyr
  laggroup2 <- rlang::ensyms(laggroup)

  #lags are computed within each laggroup level so they don't bleed
  #across group boundaries (e.g. from the end of one meeting into the
  #start of the next) -- laggroup does double duty here as both the
  #lag-reset boundary and the factor whose homogeneity is being tested
  lag.counts <- d %>%
    dplyr::select(lagcol, laggroup) %>%
    dplyr::group_by(!!!laggroup2) %>%
    dplyr::rename(lag0 = dplyr::all_of(lagcol)) %>%
    add_lags("lag0", lagnum) %>%
    dplyr::ungroup() %>%
    dplyr::group_by_all() %>%
    dplyr::summarise(freq = dplyr::n(), .groups = "drop") %>%
    as.data.frame()

  lag.counts %>%
    dplyr::select(-freq, -dplyr::all_of(laggroup)) -> lag.names

  #laggroup is treated as categorical (as.factor) since it identifies
  #discrete groups (e.g. meeting number), not a continuous predictor
  grp <- paste0("as.factor(", laggroup, ")")
  post.form <- paste(colnames(lag.names), collapse = ":")
  pre.form <- paste("freq", paste(colnames(lag.names), collapse = "*"), sep = "~")
  #saturated on the lag terms and laggroup, minus the single highest-order
  #interaction (all lags together with laggroup) -- that term has no
  #residual df to test against, same logic as shmodels()'s time12 model
  lag.form <- as.formula(paste0(pre.form, " * ", grp, " - ", post.form, ":", grp))

  options(knitr.kable.NA = '')

  lag.counts %>%
    na.omit() %>%
    glm(lag.form, family = poisson, data = .) %>%
    anova(., test = "Chisq") -> stats.out

  broom::tidy(stats.out) %>%
    knitr::kable(digits = 2, escape = F, format="html", caption = title,
                 col.names = c(
                               "Model",
                               "Model df",
                               "Deviance",
                               "Residual Deviance df",
                               "Residual Deviance",
                               "p value")) %>%
    kableExtra::kable_styling(full_width = F) %>%
    kableExtra::add_footnote(paste0("Note: Tests whether the lag structure differs across ",
                                    laggroup, " -- the highest-order lag-by-group interaction. ",
                                    "See Poole (2000)."))

  #end function
}


#internal: runs Poole's (2000, Ch. 6, "Order") order-selection procedure on
#a single ordered sequence (one group's data, or the whole sequence if
#ungrouped): "the fit of the first-order model is compared to that of the
#zeroth-order model... the fit of the second-order model [is] compared to
#that of the first-order model, the third-order model to the second-order
#model, and so on... we take the highest-order model that fits."
#
#For each depth k = 1..maxlag, a FRESHLY re-aggregated table using only
#lag0..lag_k (not the full maxlag-deep table) is fit as
#freq ~ lag0*lag1*...*lag_k, and the significance of that model's own top
#(k+1-way) interaction is tested -- this is the deviance lost by dropping
#just that term while everything else among lag0..lag_k stays in, i.e. a
#genuine test of partial association for "does order k improve on order
#k-1." Order is flagged as the highest k where this test is significant.
#
#Poole's book also describes a "quick screening device" shortcut: fit ONE
#fully saturated model on all of lag0..lag_maxlag at once, and read off
#partial-association tests for every term simultaneously (his Tables 6.5/
#6.9/6.10). That shortcut needs the saturated table to retain positive
#residual degrees of freedom to test against -- true for his large, pooled
#CIP example, but not for typical single-meeting sequences here: a 4-level
#code crossed 4 deep has far more possible cells than a single meeting has
#thought units, so the saturated model perfectly fits the *observed* cells
#regardless of which term is dropped (0 residual df for every term, so
#nothing is ever testable). The per-depth approach implemented here avoids
#that degeneracy and matches Poole's more fundamental stated logic; it also
#avoids the SB project's own Area/Area2 lesson (never leave a table split
#by lag columns the current formula doesn't use).
#noRd
order_one <- function(seq_vec, maxlag, alpha) {

  lag_names <- paste0("lag", 0:maxlag)

  build_counts <- function(k) {
    y <- data.frame(lag0 = seq_vec)
    for (i in seq_len(k)) y[[paste0("lag", i)]] <- dplyr::lag(y$lag0, i)
    y <- y[, paste0("lag", 0:k), drop = FALSE]
    y <- stats::na.omit(y)
    as.data.frame(dplyr::summarise(dplyr::group_by_all(y), freq = dplyr::n(), .groups = "drop"))
  }

  #### screen order incrementally, one depth at a time ####
  screen <- do.call(rbind, lapply(seq_len(maxlag), function(k) {
    counts.k <- build_counts(k)
    vars.k <- lag_names[1:(k + 1)]
    top.term <- paste(vars.k, collapse = ":")
    sat.form.k <- stats::as.formula(paste("freq ~", paste(vars.k, collapse = " * ")))
    sat.fit.k <- tryCatch(stats::glm(sat.form.k, data = counts.k, family = stats::poisson),
                           error = function(e) NULL)
    reduced.form.k <- stats::as.formula(paste0(
      "freq ~ ", paste(vars.k, collapse = " * "), " - ", top.term))
    reduced.fit.k <- tryCatch(stats::glm(reduced.form.k, data = counts.k, family = stats::poisson),
                               error = function(e) NULL)
    if (is.null(sat.fit.k) || is.null(reduced.fit.k)) return(NULL)
    dev.diff <- stats::deviance(reduced.fit.k) - stats::deviance(sat.fit.k)
    df.diff <- stats::df.residual(reduced.fit.k) - stats::df.residual(sat.fit.k)
    if (is.na(df.diff) || df.diff <= 0) return(NULL)
    data.frame(order = k, term = top.term, df = df.diff, partial_chisq = dev.diff,
               p_value = stats::pchisq(dev.diff, df = df.diff, lower.tail = FALSE))
  }))

  order_found <- 0L
  if (!is.null(screen)) {
    for (k in maxlag:1) {
      row <- screen[screen$order == k, ]
      if (nrow(row) == 1 && !is.na(row$p_value) && row$p_value <= alpha) {
        order_found <- k
        break
      }
    }
  }

  #confirmation that order_found is right -- both directions:
  #(a) its own term is significant (own_p, already established above)
  #(b) the NEXT order up is NOT significant, i.e. going further isn't
  #    warranted ("the third-order model did not fit" in Poole's worked
  #    example) -- already computed in screen for order_found < maxlag;
  #    at order_found == maxlag this can't be checked without raising
  #    maxlag, which is flagged rather than silently assumed
  own_row <- if (order_found > 0) screen[screen$order == order_found, ] else NULL
  own_p <- if (!is.null(own_row) && nrow(own_row) == 1) own_row$p_value else NA_real_
  next_checked <- order_found < maxlag
  next_row <- if (next_checked) screen[screen$order == order_found + 1, ] else NULL
  next_p <- if (!is.null(next_row) && nrow(next_row) == 1) next_row$p_value else NA_real_

  list(order = order_found, screen = screen, own_p = own_p,
       next_checked = next_checked, next_p = next_p)
}


#' Flag the best-fitting Markov order via tests of partial association
#'
#' Automates the order-selection procedure described in Poole (2000, Ch. 6,
#' "Order"; worked example pp. 37-38): fits the saturated model on
#' lag0..lag\code{maxlag} and computes tests of partial association (the
#' deviance increase from dropping each term while holding every other term
#' in the model, via \code{\link[stats]{drop1}}) for every term -- not the
#' sequential (Type I) tests \code{lagmodels()} reports, which can disagree
#' with partial association because the lag columns are always correlated
#' (they're lags of the same variable). The order is flagged as the highest
#' k for which the "k+1 consecutive events" term (\code{lag0:lag1:...:lag_k})
#' is still significant while the next one up is not -- "we take the
#' highest-order model that fits" (Poole, 2000, p.16). That candidate order
#' is then refit on its own (saturated on lag0..lag_k, minus its own top
#' interaction, freshly re-aggregated so unused higher lags don't silently
#' split the frequency table) and its residual-deviance goodness-of-fit is
#' reported as a confirmatory check, matching the two-step logic of Poole's
#' worked examples ("a statistical test for the fit of second-order model
#' [123] indicated good fit... the third-order model [1234] did not fit").
#'
#' This implements the two well-defined steps of Poole's procedure (partial-
#' association screening; confirmatory refit-and-test) but not the further,
#' more elaborate step in his worked examples of hand-pruning non-required
#' lower-order cross terms into a fully minimal hierarchical model -- the
#' confirmatory model here is "saturated on the needed lags minus the top
#' term," the same convention \code{lagmodels()}/\code{shmodels()}/
#' \code{hommodels()} and this project's own hand-written models already
#' use (and which has been validated against a published paper's reported
#' order/stationarity/homogeneity statistics).
#'
#' @param d dataframe containing the sequential vector
#' @param lagcol name of the column in \code{d} containing the sequential
#'   vector
#' @param laggroup name of a grouping column in \code{d} (lags reset at
#'   each level's start, and order is flagged separately per level).
#'   Default \code{""} (no grouping -- one sequence).
#' @param maxlag maximum order to screen for. Default 3.
#' @param alpha significance threshold for the partial-association tests.
#'   Default .05 -- tighten (e.g. to .01 or lower) for very large samples,
#'   where classical tests have high power to flag substantively trivial
#'   higher-order terms (Poole, 2000, pp.15-16).
#' @return Invisibly, a list (ungrouped) or data frame (grouped) with the
#'   flagged order, its own partial-association p value, and the next
#'   order's p value (confirming no need to go further); also prints a
#'   human-readable summary and the full partial-association table via
#'   \code{knitr::kable}.
#' @export

bestorder <- function(d, lagcol, laggroup = "", maxlag = 3, alpha = 0.05) {

  options(scipen = 999, warn = -1)

  if (length(d[[lagcol]]) == 0) stop("The variable does not exist in your data frame")

  `%>%` <- magrittr::`%>%`

  if (missing(laggroup) || length(laggroup) == 0 || identical(laggroup, "")) {

    res <- order_one(d[[lagcol]], maxlag = maxlag, alpha = alpha)

    cat(sprintf("Best-fitting order: %d\n", res$order))
    if (res$order > 0) {
      cat(sprintf("  order-%d term significant: p = %.4f\n", res$order, res$own_p))
    }
    if (res$next_checked) {
      cat(sprintf("  order-%d term NOT significant (confirms no need to go further): p = %.4f\n",
                  res$order + 1, res$next_p))
    } else {
      cat(sprintf("  order %d was not tested beyond maxlag = %d -- raise maxlag to check further\n",
                  res$order + 1, maxlag))
    }

    if (!is.null(res$screen)) {
      print(
        res$screen %>%
          knitr::kable(digits = 3, caption = "Tests of Partial Association (Poole, 2000)",
                       col.names = c("Order", "Term", "df", "Partial Chi-sq", "p value")) %>%
          kableExtra::kable_styling(full_width = FALSE) %>%
          kableExtra::row_spec(0, extra_css = "border-top: 1px solid black;") %>%
          kableExtra::row_spec(nrow(res$screen), extra_css = "border-bottom: 1px solid black;") %>%
          kableExtra::add_footnote(c(
            "Order: the Markov order tested (order k = does lag_k add explanatory power beyond order k-1).",
            "Term: the k+1-way interaction (lag0:lag1:...:lag_k) whose partial association is tested.",
            "df / Partial Chi-sq / p value: likelihood-ratio test of that term, holding all other order-k terms fixed (Poole, 2000).",
            sprintf("Best-fitting model: order %d (highest order whose own term is significant at alpha = %.3g).",
                    res$order, alpha)
          ), notation = "none")
      )
    }

    return(invisible(res))

  } else {

    if (length(d[[laggroup]]) == 0) stop("The group variable does not exist in your data frame")

    groups <- unique(d[[laggroup]])
    out <- lapply(groups, function(g) {
      sub <- d[d[[laggroup]] == g, , drop = FALSE]
      r <- order_one(sub[[lagcol]], maxlag = maxlag, alpha = alpha)
      data.frame(group = g, order = r$order, own_p = r$own_p,
                 next_p = r$next_p, next_checked = r$next_checked)
    })
    out <- do.call(rbind, out)
    names(out)[1] <- laggroup

    print(
      out %>%
        knitr::kable(digits = 3,
                     caption = "Best-Fitting Order by Group (Poole, 2000, partial-association procedure)",
                     col.names = c(laggroup, "Order", "Order's own p", "Next order's p",
                                   "Next order tested?")) %>%
        kableExtra::kable_styling(full_width = FALSE) %>%
        kableExtra::row_spec(0, extra_css = "border-top: 1px solid black;") %>%
        kableExtra::row_spec(nrow(out), extra_css = "border-bottom: 1px solid black;") %>%
        kableExtra::add_footnote(c(
          "Order: that group's best-fitting Markov order -- the highest order whose own defining term is significant, per Poole's (2000) tests of partial association.",
          "Order's own p: significance of the flagged order's term (order is flagged only when this is <= alpha).",
          "Next order's p: significance of the next-higher order's term -- large/non-significant confirms no need to go further; NA if not tested (see next column).",
          "Next order tested?: FALSE means Order == maxlag, so a higher order could not be checked -- raise maxlag to test further.",
          sprintf("Each row's Order value is that group's best-fitting model (alpha = %.3g).", alpha)
        ), notation = "none")
    )

    return(invisible(out))
  }

  #end function
}
