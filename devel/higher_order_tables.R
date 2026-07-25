library(lagTables)
head(lag_dat)
library(tidyverse)
summary(lag_dat)
n <- 2

lag_ll <- lag_tp <- lag_dat %>%
  rename(lag0 = seqvar)

head(lag_ll)

ipa01 <- ipa02 <- ipa_seq.dat

head(ipa01)

ipa01 <- ipa01 %>%
  rename(lag0 = area)

ipa02 <- ipa02 %>%
  rename(lag0 = area)

#https://stackoverflow.com/questions/28055927/how-can-i-automatically-create-n-lags-in-a-timeseries
#This works dynamically.  Need two versions
#Version 1 for log linear models--ascending order of lags
data.table::setDT(lag_ll)[, paste0("lag", 1:n) := data.table::shift(lag0, 1:n)][]

str(lag_ll)

#descending order of lags for transitional prob table
data.table::setDT(ipa01)[, paste0("lag", n:1) := data.table::shift(lag0, n:1)][]

#Version 1 for log linear models--ascending order of lags
data.table::setDT(ipa02)[, paste0("lag", 1:n) := data.table::shift(lag0, 1:n)][]

x <- c(rep(NA,2),1:6)
x
embed(x, 3)

#get log linear model
#how to get right variables in formula?
#Done!!!
lag_ll %>%
  as.data.frame() %>%
#  filter(study == "LAS") %>%
  dplyr::select(contains("lag")) %>%
  ftable(xtabs(, data=.)) %>%
  as.data.frame() %>%
  glm(as.formula(paste("Freq ~", paste(colnames(.)[!names(.) %in% c("Freq")],collapse="*"))), data = ., family="poisson") -> test.glm

options(warn = -1)  
anova(test.glm, test="LRT")

pchisq(deviance(test.glm), df = df.residual(test.glm), lower.tail = F)

ipa02 %>%
  as.data.frame() %>%
  filter(study == "LAS") %>%
    dplyr::select(contains("lag")) %>%
  table() %>%
  as.data.frame() %>%
  {. ->> lag.df } %>%
glm(Freq ~ lag0:lag1+lag1:lag2+lag2:lag3+lag0:lag2+lag0:lag3+lag1:lag3, data=., family = "poisson") -> test.glm1
summary(test.glm1)$deviance
anova(test.glm1, test="LRT")

#https://data.library.virginia.edu/an-introduction-to-loglinear-models/
pchisq(deviance(test.glm1), df = df.residual(test.glm1), lower.tail = F)

head(lag_dat)


n=2
las0 %>%
  rename(lag0 = area) -> las1



head(las1)

#this works well and is simple
#requires data.table, but that's okay
#data.table::setDT(las1)[, paste0("lag", 1:n) := data.table::shift(lag0, 1:n)][] %>% 
#lag_names <- paste("lag", formatC(1:1, width = nchar(max(lags)), flag = "0"), 
 #                  sep = "")

#https://gist.github.com/drsimonj/2038ff9f9c67063f384f10fac95de566
lag_functions <- setNames(paste("dplyr::lag(., ", 1:2, ")"), 
                          paste("lag", formatC(1:2, width = nchar(max(lags)), flag = "0"), 
                                sep = ""))

las1 %>% 
  mutate_at(vars(lag0), funs_(lag_functions)) %>% 
 as.data.frame() %>% 
  group_by(Condition) %>%
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


data.frame(tempcol = row.names(rbind(obs.exp, tr.std)), rbind(obs.exp, tr.std)) %>% 
    tibble::as.tibble() %>% 
    tibble::remove_rownames() %>% 
    tidyr::separate(tempcol, c("Group","Previous Unit(s)"), sep="_", extra="merge") %>%
  dplyr::mutate(`Previous Unit(s)` = str_replace(`Previous Unit(s)`, "_", "-->")) %>% 
  knitr::kable(caption=paste0("Lag ", n, "Transition Probabilities"), 
               escape = F) %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>% 
    kableExtra::pack_rows("Observed Frequencies\n(Expected Frequencies)", 1, numcodes) %>%
    kableExtra::pack_rows("Transitional Probablities\n(Standardized Residuals)", numcodes+1, numcodes*2) %>% 
  kableExtra::collapse_rows(columns = 1, valign = "top") %>% 
  kableExtra::add_header_above(c(" " = 2, "Target Unit" = ncol(obs.exp)))

  #plot
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
  
  

#add fitted values to data
cbind(lag.df,fitted(test.glm1))


