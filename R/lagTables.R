#' @importFrom stats anova as.formula chisq.test deviance df.residual ftable glm na.omit pchisq poisson xtabs
#' @importFrom utils data
NULL

#these are columns/list-elements created via dplyr/tidyr non-standard
#evaluation inside the functions below (not undefined globals) -- listed
#here so R CMD check doesn't flag them
utils::globalVariables(c(
  ".", ".data", ":=", "group", "tempcol", "Previous Unit(s)", "Var1", "Var2", "Freq",
  "freq", "data", "tidied", "time12", "model.fit", "model.stats", "fit.glanced",
  "deviance", "df.residual", "model.p", "devp.mod"
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
    expect <- DescTools::Format(as.data.frame(lag.tab$expected),leading = "drop", digits = 2)
    stdres <- DescTools::Format(as.data.frame(lag.tab$stdres), leading = "drop", digits =2)
    tr <- DescTools::Format(as.data.frame(obs/rowSums(obs)), leading = "drop", digits =2)

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
        kableExtra::pack_rows("Observed Frequencies\n(Expected Frequencies)", 1,
                              numcodes) %>%
        kableExtra::pack_rows("Transitional Probablities\n(Standardized Residuals)",
                              numcodes+1, numcodes*2) %>%
        kableExtra::collapse_rows(columns = 1, valign = "top") %>%
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
    expect <- DescTools::Format(as.data.frame(lag.tab$expected),leading = "drop", digits = 2)
    stdres <- DescTools::Format(as.data.frame(lag.tab$stdres), leading = "drop", digits =2)
    tr <- DescTools::Format(as.data.frame(obs/rowSums(obs)), leading = "drop", digits =2)

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
