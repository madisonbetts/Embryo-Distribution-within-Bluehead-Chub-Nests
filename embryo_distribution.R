#install and load all necessary packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("AICcmodavg")
#install.packages("Rtools")
install.packages("lme4")
library(readxl)
library(ggplot2)
library(dplyr)
library(AICcmodavg)
#library(Rtools)
library(lme4)

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

#print summary statistics of adjusted embyro counts by compartment 
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
a0 = glmer(cbind(sectionembryos, totnesteggs-sectionembryos) ~ (1|nest), family=binomial(link = "logit"), data=dividednests) #null model
a1 = glmer(cbind(sectionembryos, totnesteggs-sectionembryos) ~ (1|nest)+bottom+top, family=binomial(link = 'logit'), data=dividednests) #vertical gradient
a2 = glmer(cbind(sectionembryos, totnesteggs-sectionembryos) ~ (1|nest)+upstream, family=binomial(link = 'logit'), data=dividednests) #horizontal gradient
a3 = glmer(cbind(sectionembryos, totnesteggs-sectionembryos) ~ (1|nest)+MU+MD+BD+BU+TU, family=binomial(link = 'logit'), data=dividednests) #vertical and horizontal gradient

summary(a3)

cmods=list(a0,a1,a2,a3) 
Modnames=c("Null Model","Vertical Gradient Model","Horizontal Gradient Model","Vertical and Horizontal Gradient Model") 
print(aictab(cand.set = cmods, modnames = Modnames, second.ord = TRUE), digits = 2)

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
