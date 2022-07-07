#' ---
#' title: Simulations using meta-analysis
#' author: John Aponte
#' date: "2022-01-31"
#' ---

#'
#' ## Introduction
#'
#' This prigram make simulations of the trial making the analysis of each pair
#' as Agresti-Cafo and the summary using inverse meta-analysis
#'
#' ## Setup
#'
#'
#+ echo = T, result = "hide", message = F, warning = F
starttime_ <- Sys.time()
library(DBI)
library(repana)
library(plyr)
library(tidyverse)
library(readr)
library(gt)
options(knitr.duplicate.label = "allow")

#' Estimate the difference in proportions using the Agresti-Caffo CI
#'
#' @param a number of positive group 1
#' @param b total tested group 1
#' @param c number of positive group 2
#' @param d total tested group 2
#' @param id Identification for each comparison
#' @param level Level for confidence intervals
#' @return an object of class acprop with the
#' following slots:
#' p1: Proportion in group 1
#' p2: Proportion in group 2
#' dif: Differenc in proportion p2 - p1
#' var: Agresti-Caffo variance
#' lci: Lower limit for the confidence interval
#' uci: Upper limit for the confidence interval
#' w: Weight for meta-analysis
#' pw: Proportion of the weight for each observation
#' @export
acprop <- function(a,b,c,d, id, level = 0.95){
  if (missing(id)){
    id = seq(1,length(a))
  }
  p1 = a/b
  p2 = c/d
  pp1 = (a+1)/(b+2)
  pp2 = (c+1)/(d+2)
  var = pp1*(1-pp1)/(b+2) + pp2*(1-pp2)/(d+2)
  lci = (pp2-pp1) - qnorm(1-(1-level)/2)*sqrt(var)
  lci = ifelse(lci < -1, -1, lci)
  uci = (pp2-pp1) + qnorm(1-(1-level)/2)*sqrt(var)
  uci = ifelse(uci > 1, 1, uci)
  w = 1/var
  structure(
    list(
    "id" = id,
    "pos1" = a,
    "neg1" = b-a,
    "pos2" = c,
    "neg2" = d-c,
    "p1"=p1,
    "p2"=p2,
    "dif"= p2-p1,
    "var"= var,
    "lci"=lci,
    "uci"= uci,
    "w" = w,
    "pw" = w/sum(w)*100 ),
    class = "acprop"
  )
}

#' Convert to matrix an object acprop
#' @export
as.matrix.acprop <- function (x, ...){
  lx <- length(x)
  matrix(
    unlist(x[c(2:lx)]),
    ncol = lx-1,
    byrow = F,
    dimnames = list(x$id,names(x[c(2:lx)]))
  )
}

#' Convert to matrix an object acprop
#' @export
as.data.frame.acprop <- function (x, ...){
  class(x)<-"list"
  data.frame(x)
}

#' Print an acprop object
#' @export
print.acprop <- function(x, digits = 3, ...){
  print(round(as.matrix(x),digits))
  invisible(x)
}

#' Meta analysis using the inverse variance method
#'
#' @param x object with the effects and weights
#' @return and object of class metainv
#' @export
metainv <- function(x){
  UseMethod("metainv")
}

#' Default implementation for metainv
#' @export
metainv.default <- function(x, ...){
  stop(paste("Method not implemented for class", class(x)))
}

#' Implementation for metainv for objects of class acprop
metainv.acprop <- function(x, level = 0.95){
  p = sum(x$dif*x$w)/sum(x$w)
  M = length(x$w)
  # Variance according inverse weight meta-analysis
  var_m = 1/sum(x$w)
  # Variance according weighted variance method
  var_w = sum((x$dif-p)^2*x$w)/((M-1)/M*sum(x$w))
  var = var_w
  lci = p - qnorm(1-(1-level)/2)*sqrt(var)
  uci = p + qnorm(1-(1-level)/2)*sqrt(var)
  z = p/sqrt(var)
  # Homogeneity test. If significant they are not homogeneous
  q = sum(x$w*x$dif^2)- sum(x$w*x$dif)^2/sum(x$w)
  df = M-1
  # I2 measure the heterogeneity. The lower the better
  i2 = max(0,(1-(df-1)/q)*100)
  # Weighed standard deviation between cluster diff
  cv = sqrt(var)/abs(p)


  structure(
    list(
      "n" = M,
      "dif" = p,
      "var" = var,
      "lci" = lci,
      "uci" = uci,
      "z" = z,
      "prob > |z|" = (1-pnorm(abs(z)))*2,
      "Q" = q,
      "prob > Q" = pchisq(q,df,lower.tail = F),
      "I2 (%)" = i2,
      "cv" = cv
    ),
    class= "metainv"
  )
}

as.matrix.metainv <- function(x, ...){
  matrix(
    unlist(x),
    ncol = length(x),
    byrow = F,
    dimnames = list(NULL,names(x))
  )

}

ggplot_forest <- function(
  x,
  title = "Forest plot"
){
  df1 <- data.frame(as.matrix(x)) %>%
    mutate(id = 1:n()) %>%
    mutate(subjectid = row.names(.)) %>%
    mutate(invid = nrow(.)+3-id)
  df2 <- data.frame(as.matrix(metainv(x))) %>%
    mutate(invid = 1 ) %>%
    mutate(subjectid = "Combined effect")

  axisfg <-
    df1 %>% select(invid, subjectid) %>%
    bind_rows(
      df2 %>% select(invid, subjectid),
      data.frame(invid = nrow(df1)+3, subjectid = "Subject")
    ) %>%
    arrange(invid)

  nround = 2
  textdf = data.frame(
    x = 1,
    y = -1,
    textx =
      paste0(
        round(df2$dif,nround),
        " 95%CI (",
        round(df2$lci,nround),
        ", ",
        round(df2$uci,nround),
        ")"
      )
  )

  ggplot(df1) +
    aes( x= invid, y = dif, ymin = lci, ymax = uci) +
    geom_linerange( size = 0.5) +
    geom_point(aes(size= pw), shape = 15) +
    geom_linerange(data = df2, size = 0.5) +
    geom_point(data = df2, aes(size = 100), shape = 18) +
    geom_text(data = textdf,
              aes(x=x, y= y, label = textx),
              nudge_y = 0.35,
              inherit.aes = F,
              size = 0.5) +
    scale_x_continuous(breaks = axisfg$invid,labels = axisfg$subjectid) +
    scale_y_continuous("Difference in prevalence",
                       breaks = c(-1,-0.5,-0.25,0,0.25,0.5,1),
                       limits = c(-1,1),
                       expand = c(0,0)
                       )+
    coord_flip() +
    theme_minimal() +
    theme(
      axis.line = element_blank(),
      axis.line.x = element_line(size = 0.5),
      axis.ticks.x = element_line(size = 0.5),
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(15,20,5,5)
      ) +
    ggtitle(title)


}


cvia_085_optical <-
  read.csv("_data/axl05.csv") %>%
  `names<-`(., tolower(names(.)))

ooc <-
  cvia_085_optical %>%
  select(usubjid, avisit, lbanmeth, npoc,nmoc) %>%
  filter(nmoc >0) %>%
  rename(pos = npoc) %>%
  mutate(neg = nmoc-pos) %>%
  mutate(am = tolower(paste(lbanmeth,avisit, sep="_"))) %>%
  select(-nmoc, -lbanmeth, -avisit) %>%
  pivot_wider( id_cols=usubjid, names_from = "am", values_from = c("pos","neg"))



dmfa <-
  acprop(
    ooc$pos_dmfa_baseline,
    ooc$pos_dmfa_baseline+ ooc$neg_dmfa_baseline,
    ooc$pos_dmfa_final,
    ooc$pos_dmfa_final + ooc$neg_dmfa_final,
    ooc$usubjid
  )
dmfa_meta <- metainv(dmfa)

dsfa <-
  acprop(
    ooc$pos_dsfa_baseline,
    ooc$pos_dsfa_baseline+ ooc$neg_dsfa_baseline,
    ooc$pos_dsfa_final,
    ooc$pos_dsfa_final + ooc$neg_dsfa_final,
    ooc$usubjid
  )
dsfa_meta <- metainv(dsfa)

bas <-
  acprop(
    ooc$pos_dmfa_baseline,
    ooc$pos_dmfa_baseline+ ooc$neg_dmfa_baseline,
    ooc$pos_dsfa_baseline,
    ooc$pos_dsfa_baseline + ooc$neg_dsfa_baseline,
    ooc$usubjid
  )
bas_meta <- metainv(bas)

fin <-
  acprop(
    ooc$pos_dmfa_final,
    ooc$pos_dmfa_final+ ooc$neg_dmfa_final,
    ooc$pos_dsfa_final,
    ooc$pos_dsfa_final + ooc$neg_dsfa_final,
    ooc$usubjid
  )
fin_meta <- metainv(fin)

res<- rbind(
  as.matrix(dsfa_meta),
  as.matrix(dmfa_meta),
  as.matrix(bas_meta),
  as.matrix(fin_meta)
)

resdf <-
  data.frame(res) %>%
  mutate(results = c("DSFA Final vs Baseline","DMFA Final vs Baseline","Baseline DSFA vs DMFA","Final DSFA vs DMFA") ) %>%
  select(results, everything()) %>%
  rename(pvalue = prob....z.) %>%
  rename(qpvalue = prob...Q) %>%
  rename(i2 = I2....) %>%
  rename(q = Q) %>%
  mutate(pvalue = format.pval(pvalue,3,0.001)) %>%
  mutate(qpvalue = format.pval(qpvalue,3,0.001))

gt(resdf %>% select(-var,-z)) %>%
  tab_header("Table 11 Summary of combined estimations for OOcyst prevalence (Optical microscopy)(Full Analysis Population)") %>%
  fmt_number(
    columns = c("dif","lci","uci"),
    decimals = 3
  ) %>%
  fmt_number(
    columns = "i2",
    decimals = 1
  ) %>%
  fmt_number(
    columns = "cv",
    decimals = 2
  ) %>%
  cols_label(
    results = "Comparison",
    dif = "Combined\nEstimate",
    lci = "LCI",
    uci = "UCI",
    pvalue = "p_value",
    q = "Q Test",
    qpvalue = "Q Test\n p_value",
    i2 = "I2 (%)",
    cv = "CV"
  ) %>%
  cols_align("right")


ggplot_forest(dsfa, "DSFA (Final - Baseline)")
ggplot_forest(dmfa, "DMFA (Final - Baseline)")
ggplot_forest(bas, " Baseline (DSFA - DMFA)")
ggplot_forest(dmfa, "Final (DSFA - DMFA)")

#'
#' ## Analysis
#'
#+ echo = T





#'
#' ## Session Info
#'
#+  echo = F
cat("Execution date: " , format(starttime_ ),"\n")
cat("Execution time: ", round(difftime(Sys.time(),starttime_,units = "min"),2), "min\n")
print(sessionInfo(), locale = F)
