#install and load all necessary packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("AICcmodavg")
#install.packages("Rtools")
install.packages("lme4")
install.packages("MuMIn")
library(readxl)
library(ggplot2)
library(dplyr)
library(AICcmodavg)
#library(Rtools)
library(lme4)
library(MuMIn)

#read in excel sheets
dividednests <- read_excel("embryo_distribution.xlsx", sheet = "prop_eggsgrav")
View(dividednests)

genetics <- read_excel("embryo_distribution.xlsx", sheet = "genetics_analysis", na = "#N/A")
View(genetics)

#clean sheets that have NAs that interfere with analysis
summary(genetics)
genetics_cleaned <- na.omit(genetics)
View(genetics_cleaned)

options(prompt = "=")

#assign compartments as factors
dividednests$compartment <- factor(dividednests$compartment, levels=c("TD", "MD", "BD", "TU", "MU", "BU"), labels=c('1', '2', '3', '4', '5', '6'))

#let's visualize how embryo density varies across compartments
embryodensity <- ggplot(dividednests, aes(x = compartment, y = propembryoct)) +
        geom_boxplot(fill = "grey", notch = FALSE, position = position_dodge(preserve = "single")) +
        scale_x_discrete(labels = c("TD", "MD", "BD", "TU", "MU", "BU")) +  # Set x-axis labels here
        labs(
                x = "Nest Section",
                y = "Volume of Eggs (%)",
                title = ""
        ) +
        theme_bw() +
        theme(
                text = element_text(size = 15),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 15),
                plot.title = element_text(size = 15, hjust = 0.5)
        )

embryodensity

plot(dividednests$propnestvol, dividednests$propembryoct, main=" ", 
     xlab="Nest Section Volume (%)", ylab="Egg Volume (%)", pch=2, cex=1.5, cex.lab=1.5)
abline(lm(dividednests$propembryoct~dividednests$propnestvol), col="red") # regression line (y~x) 
lines(lowess(dividednests$propembryoct,dividednests$propnestvol), col="blue") # lowess line (x,y)

#let's compare sampled nest volume vs targeted volume
compare_gravvol <- ggplot(dividednests, aes(x = compartment, y = propnestvol)) +
        geom_boxplot(fill = "grey", notch = FALSE) +
        geom_hline(yintercept = 16.6, color = "blue", linewidth = 1.2, linetype = "dotted") +
        scale_x_discrete(labels = c("TD", "MD", "BD", "TU", "MU", "BU")) +
        labs(
                title = "Actual nest volumes sampled vs targeted 16.6% per section",
                x = "Nest Section",
                y = "Volume of Gravel (%)"
        ) +
        theme_bw() +
        theme(
                text = element_text(size = 15),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 15),
                plot.title = element_text(size = 15, hjust = 0.5)
        )

compare_gravvol

#plot to show how egg density varies across compartments once adjusted for sampled gravel vol
densityadjusted <- ggplot(dividednests, aes(x = compartment, y = adjembryoct)) +
        geom_boxplot(fill = "grey", notch = FALSE, position = position_dodge(preserve = "single")) +
        scale_x_discrete(labels = c("TD", "MD", "BD", "TU", "MU", "BU")) +  # x-axis labels
        labs(
                x = "Nest Section",
                y = "Adjusted Volume of Embryos (%)",
                title = ""
        ) +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                text = element_text(size = 15),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 15),
                plot.title = element_text(size = 15, hjust = 0.5)
        )

densityadjusted

#print summary statistics of adjusted embryo counts by compartment 
dividednests$compartment <- factor(dividednests$compartment,
                                   levels = c("1", "2", "3", "4", "5", "6"),
                                   labels = c("TD", "MD", "BD", "TU", "MU", "BU"))

summary_stats <- dividednests %>%
        group_by(compartment) %>%
        summarise(
                n = sum(!is.na(adjembryoct)),             # count non-missing values
                mean = mean(adjembryoct, na.rm = TRUE),
                sd = sd(adjembryoct, na.rm = TRUE),
                se = sd / sqrt(n),                        # standard error
                median = median(adjembryoct, na.rm = TRUE)
        )

print(summary_stats)

#let's investigate how embryos are distributed across nest sections, using AICcmodavg R package and adjusted embryo counts, with nest as a random effect
a0 = glmer(cbind(sectionembryos, totnestembryos-sectionembryos) ~ (1|nest), family=binomial(link = "logit"), data=dividednests) #null model
a1 = glmer(cbind(sectionembryos, totnestembryos-sectionembryos) ~ (1|nest)+bottom+top, family=binomial(link = 'logit'), data=dividednests) #vertical gradient
a2 = glmer(cbind(sectionembryos, totnestembryos-sectionembryos) ~ (1|nest)+upstream, family=binomial(link = 'logit'), data=dividednests) #horizontal gradient
a3 = glmer(cbind(sectionembryos, totnestembryos-sectionembryos) ~ (1|nest)+MU+MD+BD+BU+TU, family=binomial(link = 'logit'), data=dividednests) #vertical and horizontal gradient

summary(a3)

cmods=list(a0,a1,a2,a3) 
Modnames=c("Null Model","Vertical Gradient Model","Horizontal Gradient Model","Vertical and Horizontal Gradient Model") 
print(aictab(cand.set = cmods, modnames = Modnames, second.ord = TRUE), digits = 2)
###############################################################################
#lets remove the random effect since it is not affecting the models - this removes the issue of singularity

# Null model (no predictors)
a0_glm = glm(cbind(sectionembryos, totnestembryos - sectionembryos) ~ 1,
             family = binomial(link = "logit"),
             data = dividednests)

# Vertical gradient: bottom + top
a1_glm = glm(cbind(sectionembryos, totnestembryos - sectionembryos) ~ bottom + top,
             family = binomial(link = "logit"),
             data = dividednests)

# Horizontal gradient: upstream
a2_glm = glm(cbind(sectionembryos, totnestembryos - sectionembryos) ~ upstream,
             family = binomial(link = "logit"),
             data = dividednests)

# Vertical and horizontal gradient: MU + MD + BD + BU + TU
a3_glm = glm(cbind(sectionembryos, totnestembryos - sectionembryos) ~ MU + MD + BD + BU + TU,
             family = binomial(link = "logit"),
             data = dividednests)

summary(a3_glm)

cmods_glm=list(a0_glm,a1_glm,a2_glm,a3_glm) 
Modnames_glm=c("Null Model","Vertical Gradient Model","Horizontal Gradient Model","Vertical Horizontal Gradient Model") 
print(aictab(cand.set = cmods_glm, modnames = Modnames_glm, second.ord = TRUE), digits = 2)

library(broom)

# Now extract tidy parameter summaries with CIs for all models
model_summaries <- lapply(seq_along(cmods_glm), function(i) {
  broom::tidy(cmods_glm[[i]], conf.int = TRUE, conf.level = 0.95) %>%
    mutate(Model = Modnames_glm[i])
})

# Combine into one data frame
all_params <- dplyr::bind_rows(model_summaries) %>%
  dplyr::select(Model, term, estimate, conf.low, conf.high) %>%
  dplyr::rename(
    Variable = term,
    Estimate = estimate,
    CI_lower = conf.low,
    CI_upper = conf.high
  )

print(all_params)

# add odds ratios to the table
all_params_or <- all_params %>%
  mutate(
    OR = exp(Estimate),
    OR_lower = exp(CI_lower),
    OR_upper = exp(CI_upper)
  )

print(all_params_or)

#calculate relative variable importance as according to IT framework
cand.set <- list(
  Null_Model = a0_glm,
  Vertical_Gradient_Model = a1_glm,
  Horizontal_Gradient_Model = a2_glm,
  Vertical_Horizontal_Gradient_Model = a3_glm
)

# Calculate variable importance
sw_results <- sw(cand.set)

# View the results
print(sw_results)

#combine with parameters table
library(tibble)

# Make sure your params_tbl has a column named "Variable"
# and that sw_results is a named numeric vector

# Turn sw_results (importance scores) into a tibble
importance_tbl <- enframe(sw_results, name = "Variable", value = "Importance")

# Join the tables by Variable
final_tbl <- all_params_or %>%
  left_join(importance_tbl, by = "Variable") %>%
  relocate(Importance, .after = OR_upper)  # Put importance after ORs

# Optional: replace NA with "<0.01" for cleaner presentation
final_tbl <- final_tbl %>%
  mutate(Importance = ifelse(is.na(Importance), "<0.01", round(Importance, 2)))

# View result
print(final_tbl)
##############################################################################
#Get R² for all models using MuMIn
a_r2_mumin <- lapply(cmods, MuMIn::r.squaredGLMM)
names(a_r2_mumin) <- Modnames

for (i in seq_along(a_r2_mumin)) {
  cat(Modnames[i], ":\n")
  print(a_r2_mumin[[i]])
  cat("\n")
}

# Check for singularity
isSingular(a3)  # Should return FALSE

# Check overdispersion
a_overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, rdf = rdf, p = pval)
}
a_overdisp_fun(a3)

#let's address the issue with the random effect and overdispersion
install.packages("glmmTMB")
library(glmmTMB)

a0_tmb <- glmmTMB(cbind(sectionembryos, totnestembryos-sectionembryos) ~ (1|nest),
                  family = betabinomial(link = "logit"), data = dividednests)#null model
a1_tmb <- glmmTMB(cbind(sectionembryos, totnestembryos-sectionembryos) ~ (1|nest)+bottom+top,
                  family = betabinomial(link = "logit"), data = dividednests)#vertical gradient
a2_tmb <- glmmTMB(cbind(sectionembryos, totnestembryos-sectionembryos) ~ (1|nest)+upstream,
                  family = betabinomial(link = "logit"), data = dividednests)#horizontal gradient
a3_tmb <- glmmTMB(cbind(sectionembryos, totnestembryos - sectionembryos) ~ MU + MD + BD + BU + TU + (1|nest),
                  family = betabinomial(link = "logit"), data = dividednests)#vertical and horizontal gradient

summary(a3_tmb)
cmods=list(a0_tmb,a1_tmb,a2_tmb,a3_tmb) 
Modnames_tmb=c("Null Model","Vertical Gradient Model","Horizontal Gradient Model","Vertical and Horizontal Gradient Model") 
print(aictab(cand.set = cmods, modnames = Modnames_tmb, second.ord = TRUE), digits = 2)

# Marginal/conditional R² for glmmTMB
library(performance)
r2(a3_tmb)
#Check fixed-effect collinearity
library(performance)
check_collinearity(a3_tmb) #want VIF<5

#inspect fixed effects for NA estimates
summary(a3_tmb)$varcor

# Check overdispersion
a_tmb_overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, rdf = rdf, p = pval)
}
a_overdisp_fun(a3)

# Confidence intervals
confint(a3_tmb)
################################################################################
#now let's look into how host, associate, and manipulator embryos are distributed across a nest

#plot how density of host embryos changes by nest section
hosteggs <- ggplot(genetics_cleaned, aes(x=compartment, y=host_perc)) + geom_boxplot(position = position_dodge(preserve = "single"))
hosteggs + labs(x="Nest Section", y="Host Eggs (%)", title="Proportion of host eggs by nest compartment")

hostembryo <- ggplot(genetics_cleaned, aes(x = compartment, y = host_perc)) +
        geom_boxplot(fill = "grey", notch = FALSE, position = position_dodge(preserve = "single")) +
        scale_x_discrete(labels = c("TD", "MD", "BD", "TU", "MU", "BU")) +
        labs(
                x = "Nest Compartment",
                y = "Host Embryo (%)",
                title = "Proportion of host embryo by nest compartment"
        ) +
        theme_bw() +
        theme(
                text = element_text(size = 15),
                axis.title = element_text(size = 15),
                axis.text = element_text(size = 15),
                plot.title = element_text(size = 15, hjust = 0.5)
        )

hostembryo

#plot how density of manipulator embryos change by nest compartment
manipembryo <- ggplot(genetics_cleaned, aes(x=compartment, y=manip_perc)) + 
        geom_boxplot(fill = "grey", notch = FALSE, position = position_dodge(preserve = "single")) +
  scale_x_discrete(labels = c("TD", "MD", "BD", "TU", "MU", "BU")) +
  labs(
    x = "Nest Compartment",
    y = "Manipulator Embryo (%)",
    title = "Proportion of manipulator embryo by nest compartment"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    plot.title = element_text(size = 15, hjust = 0.5)
  )

manipembryo


#let's investigate how host/associate and manipulator/nonmanipulator embryos are distributed across nest compartments
b0 = glmer(cbind(host_total, associate_total) ~ (1|nest), family=binomial(link='logit'), data=genetics_cleaned)
b1 = glmer(cbind(host_total, associate_total) ~ (1|nest)+MD+MU+BD+BU+TU, family=binomial(link='logit'), data=genetics_cleaned)
b2 = glmer(cbind(host_total, associate_total) ~ (1|nest)+bottom+top, family=binomial(link='logit'), data=genetics_cleaned)
b3 = glmer(cbind(host_total, associate_total) ~ (1|nest)+upstream, family=binomial(link='logit'), data=genetics_cleaned)

c0 = glmer(cbind(manipulator_total, nonmanipulator_total) ~ (1|nest), family=binomial(link='logit'), data=genetics_cleaned)
c1 = glmer(cbind(manipulator_total, nonmanipulator_total) ~ (1|nest)+MD+MU+BD+BU+TU, family=binomial(link='logit'), data=genetics_cleaned)
c2 = glmer(cbind(manipulator_total, nonmanipulator_total) ~ (1|nest)+bottom+top, family=binomial(link='logit'), data=genetics_cleaned)
c3 = glmer(cbind(manipulator_total, nonmanipulator_total) ~ (1|nest)+upstream, family=binomial(link='logit'), data=genetics_cleaned)

cmods2=list(b0,b1,b2,b3)
Modnames=c("Null Model","Vertical Gradient Model","Horizontal Gradient Model","Vertical and Horizontal Gradient Model") 
print(aictab(cand.set = cmods2, modnames = Modnames, second.ord = TRUE), digits = 2)


cmods3=list(c0,c1,c2,c3) 
Modnames=c("Null Model","Vertical Gradient Model","Horizontal Gradient Model","Vertical and Horizontal Gradient Model") 
print(aictab(cand.set = cmods3, modnames = Modnames, second.ord = TRUE), digits = 2)

####### Larvae vs Embryo Analysis #######
#read in larvae data
larv <- read_excel("embryo_distribution.xlsx", sheet = "larv_analysis")
View(larv)

#install packages for larvae analysis
install.packages("lmerTest")
install.packages("emmeans")
install.packages("pbkrtest")
library(lmerTest)
library(emmeans)
library(pbkrtest)

#let's compare larvae and egg proportions across nest sections using a linear mixed effects model with nest as a random effect
larvmodel <- lmer(larv_prop ~ compartment + (1 | nest), data = larv)
summary(larvmodel)

#produce pairwise comparisons for all compartments using Tukey HSD post-hoc
larv_pairwise <- emmeans(larvmodel, pairwise ~ compartment, adjust = "tukey")
larv_pairwise

#produce table with upper and lower bounds of confidence interval for all pairwise comparisons
larv_comparisons <- larv_pairwise$contrasts
larv_comparisions_summary <- summary(larv_pairwise$contrasts, infer = c(TRUE, TRUE))
larv_comparisions_summary[, c("contrast", "estimate", "lower.CL", "upper.CL")]
