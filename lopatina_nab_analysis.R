library(tidyverse)
library(ggpubr)
library(rstatix)
library(plyr)

# data preprocessing:
# turn numerical variable into factors, if needed, add columns with log(Nab) and log(IgG)
# add factor variable, which contains two levels of donors' Nab concentration
# merge recipients and donors data into one dataframe


donors <- read.csv(file = "Nab_donors_IB.csv",
                   header = TRUE)
donors <- data.frame(donors)

for (var in 1:length(donors$Nab_donor)) {
  if (donors$Nab_donor[var] == 0) {
    donors$Nab_donorlg[var] <- 0
    donors$Nab_donor[var] = 1 
  } else {
    donors$Nab_donorlg[var] <- log(donors$Nab_donor[var], 10)
  }
}


donors$IgG_donor <- as.numeric(donors$IgG_donor)
donors$ID_donor <- as.factor(donors$ID_donor)
donors$ID_recipient <- as.factor(donors$ID_recipient)

for (var in 1:length(donors$IgG_donor)) {
  if (donors$IgG_donor[var] == 0) {
    donors$IgG_donorlg[var] <- 0
  } else {
    donors$IgG_donorlg[var] <- log(as.numeric(donors$IgG_donor[var]), 10)
  }
}

recipients <- read.csv(file = "Nab_recipients_IB.csv",
                       header = TRUE)
recipients <- data.frame(recipients)

recipients$sampling_spot <- as.factor(recipients$sampling_spot)
recipients$ID <- as.factor(recipients$ID)
recipients$IgG <- as.numeric(recipients$IgG)

for (var in 1:length(recipients$Nab)) {
  if (recipients$Nab[var] == 0) {
    recipients$Nab_lg[var] <- 0
    recipients$Nab[var] = 1
  } else {
    recipients$Nab_lg[var] <- log(recipients$Nab[var],10)
  }
}

for (var in 1:length(recipients$IgG)) {
  if (recipients$IgG[var] == 0) {
    recipients$IgG_lg[var] <- 0
    recipients$IgG[var] = 1
  } else {
    recipients$IgG_lg[var] <- log(recipients$IgG[var],10)
  }
}

joined <- merge(x=recipients, y=donors, by.x="ID", by.y="ID_recipient")

for (var in 1:length(joined$Nab_donor)) {
  if (joined$Nab_donor[var] < 80) {
    joined$donors_level[var] <- "<80"
  } else {
    joined$donors_level[var] <- ">=80"
  }
}

joined$donors_level <- as.factor(joined$donors_level)

# get group charactestic

summary(joined)


# general analysis of Nab increase (indepent variable: sampling_spot (time)) via anova and paired wilcoxon test

res <- wilcox.test(joined$sampling_spot, joined$Nab_lg, paired = TRUE)
res

aov_ <- aov(Nab ~ donors_level*sampling_spot + Error(ID/sampling_spot), data=joined)
summary(aov_)
model.tables(aov_, "means")

compare_means(Nab_lg ~ sampling_spot,  data = joined, method = "anova")

ggboxplot(joined, x = "sampling_spot", y = "Nab_lg",
          fill = "sampling_spot", palette = "jco")+
  stat_compare_means()+theme_bw()+ggtitle("Нарастание нейтрализующих антител у реципиентов после трансфузии")


ggboxplot(joined, x = "sampling_spot", y = "Nab_lg",
          fill = "sampling_spot", palette = "jco")+
  stat_compare_means(method = "anova")+theme_bw()+ggtitle("Нарастание нейтрализующих антител у реципиентов после трансфузии")

# pairwise comparisons (using sampling_spot == 0 group as reference)

compare_means(Nab_lg ~ sampling_spot,  data = joined)
my_comparisons <- list( c("0", "1"), c("0", "2"), c("0", "3"), c("0", "4"))


ggboxplot(joined, x = "sampling_spot", y = "Nab_lg",
          fill = "sampling_spot", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means()     # Add global p-value


ggboxplot(joined, x = "sampling_spot", y = "Nab_lg",
          fill = "sampling_spot", palette = "jco")+
  stat_compare_means(method = "anova", label.y = 3)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "0") 



# Analysis of Nab increase by two variables (sampling spot and donors' Nab concentration); p-value via wilcoxon test included

ggboxplot(data = joined, x = "donors_level", y = "Nab", fill = "donors_level", title = "Зависимость титра нейтрализующих антител от момента забора и насыщенности плазмы", xlab = "Насыщенность плазмы доноров", ylab ="Уровень нейтрализующих антител реципиентов")+facet_wrap(joined$sampling_spot)+scale_y_log10()+theme_bw()+stat_compare_means()


# Two-way mixed RM-ANOVA

joined %>%
  group_by(sampling_spot, donors_level) %>%
  get_summary_stats(Nab_lg, type = "mean_sd")

bxp <- ggboxplot(
  joined, x = "sampling_spot", y = "Nab",
  fill = "donors_level")
bxp

joined %>%
  group_by(sampling_spot, donors_level) %>%
  identify_outliers(Nab)

joined %>%
  group_by(sampling_spot, donors_level) %>%
  shapiro_test(Nab_lg)

ggqqplot(joined, "Nab_lg", ggtheme = theme_bw()) +
  facet_grid(sampling_spot ~ donors_level)

joined %>%
  group_by(sampling_spot) %>%
  levene_test(Nab_lg ~ donors_level)

box_m(joined[, "Nab_lg", drop = FALSE], joined$donors_level)

res.aov <- anova_test(
  data = joined, dv = Nab_lg, wid = ID,
  between = donors_level, within = sampling_spot
)
get_anova_table(res.aov)


one.way <- joined %>%
  group_by(sampling_spot) %>%
  anova_test(dv = Nab, wid = ID, between = donors_level) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


pwc <- joined %>%
  group_by(sampling_spot) %>%
  pairwise_t_test(Nab ~ donors_level, p.adjust.method = "bonferroni")
pwc

#

ggboxplot(data = joined, x = "donors_level", y = "Nab", fill = "donors_level", title = "Зависимость титра нейтрализующих антител от момента забора и насыщенности плазмы", xlab = "Насыщенность плазмы доноров", ylab ="Уровень нейтрализующих антител реципиентов")+facet_wrap(joined$sampling_spot)+scale_y_log10()+theme_bw()+stat_compare_means()
ggboxplot(data = joined, x = "sampling_spot", y = "Nab", fill = "donors_level", title = "Зависимость титра нейтрализующих антител от момента забора и насыщенности плазмы", xlab = "Момент времени после трансфузии", ylab ="Уровень нейтрализующих антител реципиентов")+scale_y_log10()+theme_bw()+stat_compare_means(method = "wilcox.test", paired = T)

#debugs for MANOVA 

debug_contr_error <- function (dat, subset_vec = NULL) {
  if (!is.null(subset_vec)) {
    ## step 0
    if (mode(subset_vec) == "logical") {
      if (length(subset_vec) != nrow(dat)) {
        stop("'logical' `subset_vec` provided but length does not match `nrow(dat)`")
      }
      subset_log_vec <- subset_vec
    } else if (mode(subset_vec) == "numeric") {
      ## check range
      ran <- range(subset_vec)
      if (ran[1] < 1 || ran[2] > nrow(dat)) {
        stop("'numeric' `subset_vec` provided but values are out of bound")
      } else {
        subset_log_vec <- logical(nrow(dat))
        subset_log_vec[as.integer(subset_vec)] <- TRUE
      } 
    } else {
      stop("`subset_vec` must be either 'logical' or 'numeric'")
    }
    dat <- base::subset(dat, subset = subset_log_vec)
  } else {
    ## step 1
    dat <- stats::na.omit(dat)
  }
  if (nrow(dat) == 0L) warning("no complete cases")
  ## step 2
  var_mode <- sapply(dat, mode)
  if (any(var_mode %in% c("complex", "raw"))) stop("complex or raw not allowed!")
  var_class <- sapply(dat, class)
  if (any(var_mode[var_class == "AsIs"] %in% c("logical", "character"))) {
    stop("matrix variables with 'AsIs' class must be 'numeric'")
  }
  ind1 <- which(var_mode %in% c("logical", "character"))
  dat[ind1] <- lapply(dat[ind1], as.factor)
  ## step 3
  fctr <- which(sapply(dat, is.factor))
  if (length(fctr) == 0L) warning("no factor variables to summary")
  ind2 <- if (length(ind1) > 0L) fctr[-ind1] else fctr
  dat[ind2] <- lapply(dat[ind2], base::droplevels.factor)
  ## step 4
  lev <- lapply(dat[fctr], base::levels.default)
  nl <- lengths(lev)
  ## return
  list(nlevels = nl, levels = lev)
}

debug_contr_error2 <- function (form, dat, subset_vec = NULL) {
  ## step 0
  if (!is.null(subset_vec)) {
    if (mode(subset_vec) == "logical") {
      if (length(subset_vec) != nrow(dat)) {
        stop("'logical' `subset_vec` provided but length does not match `nrow(dat)`")
      }
      subset_log_vec <- subset_vec
    } else if (mode(subset_vec) == "numeric") {
      ## check range
      ran <- range(subset_vec)
      if (ran[1] < 1 || ran[2] > nrow(dat)) {
        stop("'numeric' `subset_vec` provided but values are out of bound")
      } else {
        subset_log_vec <- logical(nrow(dat))
        subset_log_vec[as.integer(subset_vec)] <- TRUE
      } 
    } else {
      stop("`subset_vec` must be either 'logical' or 'numeric'")
    }
    dat <- base::subset(dat, subset = subset_log_vec)
  }
  ## step 0 and 1
  dat_internal <- stats::lm(form, data = dat, method = "model.frame")
  attr(dat_internal, "terms") <- NULL
  ## rely on `debug_contr_error` for steps 2 to 4
  c(list(mf = dat_internal), debug_contr_error(dat_internal, NULL))
}

NA_preproc <- function (dat) {
  for (j in 1:ncol(dat)) {
    x <- dat[[j]]
    if (is.factor(x) && anyNA(x)) dat[[j]] <- base::addNA(x)
    if (is.character(x)) dat[[j]] <- factor(x, exclude = NULL)
  }
  dat
}

info <- debug_contr_error2(Nab_lg ~ ., joined)
dat <- info$mf
nrow(dat)
info$nlevels

new_joined <- NA_preproc(joined)

new_info <- debug_contr_error2(Nab_lg ~ ., new_joined)

new_info$nlevels

# new life


ggplot(data = joined, aes(x = sampling_spot, y = Nab, group = sampling_spot))+geom_boxplot(fill="gray")+scale_y_log10()+theme_bw()+ggtitle("Нарастание нейтрализующих антител у реципиентов после трансфузии")+
  theme(plot.title = element_text(hjust = 0.5))



