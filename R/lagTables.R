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
    
    lagdat <- d %>%
      dplyr::rename(lag0 = lagvar, group=laggroup)
    
    #https://gist.github.com/drsimonj/2038ff9f9c67063f384f10fac95de566
    #gets lags without using data.table...
    lag_functions <- setNames(paste("dplyr::lag(., ", 1:lagnum, ")"), 
                              paste("lag", formatC(1:lagnum, 
                                                   width = nchar(max(1:lagnum)), 
                                                   flag = "0"), sep = ""))
    lagdat %>% 
      mutate_at(vars(lag0), funs_(lag_functions)) %>% 
      as.data.frame() %>% 
      group_by(group) %>%
      dplyr::select(contains("lag")) %>%
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
        tibble::as.tibble() %>% 
        tibble::remove_rownames() %>% 
        tidyr::separate(tempcol, c("Group","Previous Unit(s)"), sep="_", extra="merge") %>%
        dplyr::mutate(`Previous Unit(s)` = str_replace(`Previous Unit(s)`, "_", "-->")) %>% 
        knitr::kable(caption=paste0("Lag ", n, "Transition Probabilities By Group"), 
                     escape = F) %>%
        kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% 
        kableExtra::pack_rows("Observed Frequencies\n(Expected Frequencies)", 1, numcodes) %>%
        kableExtra::pack_rows("Transitional Probablities\n(Standardized Residuals)", numcodes+1, numcodes*2) %>% 
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
          group_by(Var1) %>% 
          ggplot2::ggplot(ggplot2::aes(x=Var2, y=Freq, fill=Var2)) +
          ggplot2::geom_bar(stat='identity', width = .1) +
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
      dplyr::rename(lag0 = lagvar)

    #https://gist.github.com/drsimonj/2038ff9f9c67063f384f10fac95de566
    lag_functions <- setNames(paste("dplyr::lag(., ", 1:lagnum, ")"), 
                              paste("lag", formatC(1:lagnum, 
                                                   width = nchar(max(1:lagnum)), 
                                                   flag = "0"), sep = ""))
    lagdat %>% 
      mutate_at(vars(lag0), funs_(lag_functions)) %>% 
      as.data.frame() %>% 
      #group_by(Condition) %>%
      dplyr::select(contains("lag")) %>%
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
    
    #data.frame(tempcol = row.names(obs.exp), obs.exp) %>% 
     # separate(tempcol, c("Group","Comment"))
    
    numcodes <- nrow(obs.exp)
    
    print(
      data.frame(tempcol = row.names(rbind(obs.exp, tr.std)), rbind(obs.exp, tr.std)) %>% 
        tibble::as.tibble() %>% 
        tibble::remove_rownames() %>%
        dplyr::mutate(tempcol = str_replace(tempcol, "_", "-->")) %>% 
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
  
  
  #create data frame with lags = lagnum
  #https://gist.github.com/drsimonj/2038ff9f9c67063f384f10fac95de566
  lags <- seq(lagnum)
  lag_names <- paste("lag", formatC(lags, width = nchar(max(lags)), flag = "0"), 
                     sep = "")
  lag_functions <- setNames(paste("dplyr::lag(., ", lags, ")"), lag_names)
  
  #define magrittr pipe
  `%>%` <- magrittr::`%>%`
  
  if(missing(laggroup)) {
    
    #print("In non-group loop")
    
    lag.counts <- d %>%
      dplyr::select(lagcol) %>%
      dplyr::mutate_at(dplyr::vars(lagcol), dplyr::funs_(lag_functions)) %>%
      dplyr::rename(lag0 = lagcol) %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = n()) %>%
      dplyr::ungroup() %>%
      as.data.frame()
    
    #the ncol routine removes freq from the end of the formula
    lag.form <- as.formula(paste("freq", paste(colnames(lag.counts[-(ncol(lag.counts))]), 
                                               collapse=" * "), sep="~"))
    
    fit <- glm(lag.form, data=lag.counts, family=poisson)
    
    fit.out <- anova(fit, test="Chisq")
    
    #prints to console if requested
    #suppressWarnings(print(broom::tidy(fit.out)))
    
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
    
    #Use this loop if group condition specified
  } else {
    
    if(length(laggroup)==0) {stop("The group variable does not exist in your data frame")}
    
    lag.counts <- d %>%
      dplyr::select(lagcol, laggroup) %>%
      dplyr::mutate_at(dplyr::vars(lagcol), dplyr::funs_(lag_functions)) %>%
      dplyr::rename(lag0 = lagcol) %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = n()) %>%
      dplyr::ungroup() %>%
      as.data.frame()
    
    #print(length(laggroup))
    
    lag.counts %>%
      dplyr::select(-freq, -laggroup) -> lag.names
    
    lag.form <- as.formula(paste("freq", paste(colnames(lag.names), collapse=" * "), sep="~"))
    
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
      knitr::kable(digits = 2, 
                   col.names = c(laggroup,
                                 "Model",
                                 "Model df",
                                 "Deviance",
                                 "Residual Deviance df",
                                 "Residual Deviance",
                                 "p value")) %>%
      kableExtra::kable_styling(full_width = F)
      #kableExtra::collapse_rows(columns = 1:length(laggroup), valign = "top")
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

shmodels <- function(d, lagcol, laggroup="", title="Homogeniety Tests for Log Linear Models", lagnum=1, st=1) {

  #https://data.library.virginia.edu/an-introduction-to-loglinear-models/
  
  options(scipen = 999, warn = -1)
  
  #make sure sequential data column exists
  if(length(d[[lagcol]])==0) {stop("The variable does not exist in your data frame")}
  
  #create data frame with lags = lagnum
  #https://gist.github.com/drsimonj/2038ff9f9c67063f384f10fac95de566
  lags <- seq(lagnum)
  lag_names <- paste("lag", formatC(lags, width = nchar(max(lags)), flag = "0"), 
                     sep = "")
  lag_functions <- setNames(paste("dplyr::lag(., ", lags, ")"), lag_names)
  
  #define magrittr pipe
  `%>%` <- magrittr::`%>%`
  
  if(missing(laggroup)) {
  
    print("You made it this far")
    lag.counts <- d %>%
      dplyr::select(lagcol) %>%
      dplyr::mutate(time12 = ifelse(row_number() < n()/2, 1,2)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate_at(dplyr::vars(lagcol), dplyr::funs_(lag_functions)) %>%
      dplyr::rename(lag0 = lagcol) %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = n()) %>%
      dplyr::ungroup() %>%
      as.data.frame()
    
    lag.counts %>%
      dplyr::select(-freq, -time12) -> lag.names
    
    post.form <- paste(colnames(lag.names), collapse  = ":")
    print(post.form)
    
    post0.form <- paste(colnames(lag.names[-length(lag.names)]), collapse  = ":")
    print(post0.form)
    
    pre.form <- paste("freq", paste(colnames(lag.names), collapse="*"), sep="~")
    print(pre.form)
    
    lag.form <- as.formula(paste0(pre.form, "*time12", " - ",post.form,":time12"))
    print(lag.form)
    
    ant.form <- paste0("freq", "~", post.form, " + ", post0.form, ":time12")
    
    print(ant.form)
    
    model.form <- as.formula(ant.form)
    print(model.form)

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
    
    } else {
  
      #Make sure grouping variable exists  
  if(length(d[[laggroup]])==0) {stop("The group variable does not exist in your data frame")}
      
  #handles a weird quoting issue with multiple calls to
  #the same var in dplyr
  laggroup2 <- ensyms(laggroup)

    lag.counts <- d %>%
      dplyr::select(lagcol, laggroup) %>%
      #add time to df if stat model requested
      dplyr::group_by(!!!laggroup2) %>%
      dplyr::mutate(time12 = ifelse(row_number() < n()/2, 1,2)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate_at(dplyr::vars(lagcol), dplyr::funs_(lag_functions)) %>%
      dplyr::rename(lag0 = lagcol) %>%
      dplyr::group_by_all() %>%
      dplyr::summarise(freq = n()) %>%
      dplyr::ungroup() %>%
      as.data.frame()
      
    #print(length(laggroup))
     
    lag.counts %>%
      dplyr::select(-freq, -laggroup, -time12) -> lag.names
    
    post.form <- paste(colnames(lag.names), collapse  = ":")
    print(post.form)
    
    post0.form <- paste(colnames(lag.names[-length(lag.names)]), collapse  = ":")
    print(post0.form)
    
    pre.form <- paste("freq", paste(colnames(lag.names), collapse="*"), sep="~")
    print(pre.form)
    
    lag.form <- as.formula(paste0(pre.form, "*time12", " - ",post.form,":time12"))
    print(lag.form)
    
    ant.form <- paste0("freq", "~", post.form, " + ", post0.form, ":time12")
    
    print(ant.form)
    
    model.form <- as.formula(ant.form)
    print(model.form)
    
    
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

    test.out %>%
      tidyr::unnest(model.stats, .drop=T) -> stats.out
    
    test.out %>%
     tidyr::unnest(fit.glanced, .drop=T) %>%
      dplyr::mutate(model.p = pchisq(deviance, df = df.residual, lower.tail = F)) %>%
      dplyr::mutate(devp.mod = paste0(round(deviance,2),"<br>(p. = ",round(model.p),2,")")) %>%
      dplyr::select(laggroup, devp.mod) -> fit.out

    #print(dplyr::left_join(fit.out, stats.out))
    #stop()
    
        
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

