#' Create Transitional Probability Tables and Plots
#' 
#' This function works with a vector of sequential data and creates
#' html tables for the observed and expected frequencies, as well
#' as the transitional probabilities and standardized residuals
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

trprobs <- function(d, lagvar, laggroup=NULL, title="Lag Sequential Descriptive Stats", lagnum=1, plots=0, dname="") {
  
  options(scipen = 999)
  
  d <- as.data.frame(d, stringsAsFactors=FALSE)
  
  if(length(d[lagvar])==0) {stop("The variable does not exist in your data frame")}
  
  `%>%` <- magrittr::`%>%`
  
  #get number of unique codes in dataset
  d %>%
    dplyr::select(lagvar) %>%
    na.omit() %>%
    dplyr::n_distinct(.) -> numcodes
  
  d %>%
    na.omit() %>%
    dplyr::arrange_at(., lagvar) %>%
    dplyr::select(lagvar) %>%
    dplyr::distinct() %>%
    unlist() -> codenames
  
  #get better facet labels for plots
  better_label <- function(string) {
    string <- paste0("Comments following ",string)
  }
  
  #IF GROUPING VARIABLE
  if(length(d[laggroup])!=0) {
    
    #get number of groups
    d %>%
      dplyr::select(laggroup) %>%
      na.omit() %>%
      dplyr::n_distinct(.) -> numgroupcodes
    
    #get group names
    d %>%
      na.omit() %>%
      dplyr::arrange_at(., laggroup) %>%
      dplyr::select(laggroup) %>%
      dplyr::distinct() %>%
      unlist() -> groupnames
    
    #gets dynamic 2nd header info
    #see last comment in: https://stackoverflow.com/questions/45206908/kableextra-dynamic-add-header-above-labeling
    groupnames <- sub("(.)", "\\U\\1", groupnames, perl=TRUE)
    tablecols <- rep(numcodes,2)
    names(tablecols) <- c(groupnames)
    
    lagdat <- d %>%
      dplyr::select(lagvar, laggroup) %>%
      dplyr::rename(raw = lagvar) %>%
      dplyr::mutate(lag = dplyr::lag(raw, lagnum)) 
    
    print(
      lagdat %>%
        dplyr::select(lag, raw, laggroup) %>%
        table() %>%
        apply(3, chisq.test, simulate.p.value = TRUE) %>%
        lapply(`[`, c(9)) %>%
        reshape2::melt() %>%
        tidyr::spread(key = L2, value = value) %>%
        dplyr::rename(Condition = L1) %>%
        dplyr::select(Condition, lag, raw, stdres) %>%
        dplyr::arrange(Condition) %>%
        ggplot2::ggplot(aes(x=lag, y=stdres, fill=Condition)) +
        ggplot2::geom_bar(stat='identity', position=position_dodge(), width = .5) +
        ggplot2::geom_hline(yintercept = 1.96, linetype="dashed") +
        ggplot2::geom_hline(yintercept = -1.96, linetype="dashed") +
        ggplot2::coord_flip() +
        ggplot2::facet_wrap(~ raw) +
        ggplot2::scale_fill_grey() +
        ggplot2::labs(y="Standardized Residuals", x = "", fill = "Condition",
                title = paste0("Standardized Residuals For All Comment Types, Data = ", 
                               dname, ", Lag = ",lagnum)))
    
    #print table to viewer
    lagdat %>% 
      dplyr::mutate(lag1 = lag(lagvar)) %>%
      dplyr::select(lagvar, lag1, Condition) %>%
      table() %>%
      apply(3, chisq.test, simulate.p.value = TRUE) %>%
      lapply(`[`, c(6,7,9)) %>%
      reshape2::melt() %>%
      tidyr::spread(key = L2, value = value) %>%
      dplyr::rename(Condition = L1) %>%
      dplyr::arrange(Condition, lagvar, lag1) %>%
      dplyr::group_by(Condition, lagvar) %>%
      dplyr::mutate(trprob = observed/sum(observed)) %>%
      dplyr::mutate(obsexp = paste0(observed,"<br>(",round(expected,2),")")) %>%
      dplyr::mutate(tpsres = paste0(round(stdres,2),"<br>(",round(trprob,2),")")) %>%
      assign("stats.out",.,envir = .GlobalEnv) %>%
      dplyr::select(lag1, lagvar, Condition, obsexp) %>%
      reshape2::acast(., lagvar ~ lag1 ~ Condition ) %>%
      tbl_df() %>% 
      as.data.frame() -> top.dat
    
    stats.out %>%
      dplyr::select(lag1, lagvar, Condition, tpsres) %>%
      reshape2::acast(., lagvar ~ lag1 ~ Condition ) %>%
      tbl_df() %>% 
      as.data.frame() -> bot.dat
    
    lagout <- rbind(top.dat, bot.dat)
    lagout <- cbind(rep(codenames, 2),lagout)
    
    lagout %>%
      knitr::kable(caption = paste0(title, ", Lag = ", lagnum), digits=2, format = "html", escape = F, align=rep('c', numcodes), col.names = c("",rep(codenames,2))) %>%
      kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE)  %>%
      kableExtra::add_header_above(c("", tablecols)) %>%
      kableExtra::group_rows("Observed Frequencies\n(Expected Frequencies)", 1, numcodes) %>%
      kableExtra::group_rows("Transitional Probablities\n(Standardized Residuals)", numcodes+1, numcodes*2)
    
    #IF NO GROUP VARIABLE
  } else {
    
    lagdat <- d %>%
      dplyr::select(lagvar) %>%
      dplyr::rename(raw = lagvar) %>%
      dplyr::mutate(lag = dplyr::lag(raw, lagnum)) 
    
    #function for trans prob matrixy
    trans.matrix <- function(X, prob=T)
    {
      tt <- X
      if(prob) tt <- tt / rowSums(tt)
      round(tt,2)
    }
    
    transp <- trans.matrix(chisq.test(table(lagdat$lag, lagdat$raw), simulate.p.value = TRUE)$observed)
    obs <- chisq.test(table(lagdat$lag, lagdat$raw), simulate.p.value = TRUE)$observed
    exp <- chisq.test(table(lagdat$lag, lagdat$raw), simulate.p.value = TRUE)$expected
    sres <- chisq.test(table(lagdat$lag, lagdat$raw), simulate.p.value = TRUE)$stdres
    
    obexp <- rbind(as.data.frame.matrix(obs), 
                   DescTools::Format(as.data.frame.matrix(transp),
                                     leading = "drop", digits = 2))
    
    psres <- rbind(DescTools::Format(as.data.frame.matrix(exp), digits=2), 
                   DescTools::Format(as.data.frame.matrix(sres), digits=2))
    
    mypaste <- function(x,y) paste0(x, "<br>(", y, ")")
    
    lagout <- mapply(mypaste, obexp, psres)
    
    rownames(lagout) <- rep(colnames(obs),2)
    
    #provide plots if requested
    if(plots>0) {  
      print(
        sres %>%
          as.data.frame() %>%
          dplyr::group_by(Var1) %>%
          ggplot2::ggplot(aes(x=Var2, y=Freq, fill=Var2)) +
          ggplot2::geom_bar(stat='identity', width = .5) +
          ggplot2::geom_hline(yintercept = 1.96, linetype="dashed") +
          ggplot2::geom_hline(yintercept = -1.96, linetype="dashed") +
          ggplot2::coord_flip() +
          ggplot2::facet_wrap(~ Var1, labeller = labeller(Var1=better_label)) +
          ggplot2::scale_fill_grey() +
          ggplot2::labs(y="Standardized Residuals", x = "", fill = "Comment\nType",
                        title = paste0("Standardized Residuals For All Comment Types, Data = ", dname, ", Lag = ",lagnum))
      )
    }
    
    #print table to viewer
    lagout %>% 
      knitr::kable(caption = paste0(title, ", Lag = ", lagnum), digits=2, format = "html", escape = F, align=rep('c', numcodes), col.names = c(rep(codenames,1))) %>%
      kableExtra::kable_styling(bootstrap_options = "striped", full_width = FALSE)  %>%
      kableExtra::group_rows("Observed Frequencies\n(Expected Frequencies)", 1, numcodes) %>%
      kableExtra::group_rows("Transitional Probablities\n(Standardized Residuals)", numcodes+1, numcodes*2)
  }
  #End function 
}


#' Estimates log linear models for sequential data
#' Requires a vector of sequential data
#' Can request estimates by group
#' @param d dataframe containing sequential vector
#' @param lagcol vector in dataframe containing sequential vector
#' @param laggroup vector in dataframe containing grouping variable.  Default is empty (no groups)
#' @param lagnum number of lags. Default is 1
#' @param title caption for table
#' @return Table with log linear estimates for lag sequential data
#' @export

lagmodels <- function(d, lagcol, laggroup, title="Lag Sequential Log Linear Models", lagnum=1) {
  
  options(scipen = 999)
  
  if(length(lagcol)==0) {stop("The variable does not exist in your data frame")}
  
  
  #create data frame wiht lags = lagnum
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
      kableExtra::kable_styling(full_width = F) %>%
      kableExtra::collapse_rows(columns = 1:length(laggroup), valign = "top")
  }
  
  #end function
}
