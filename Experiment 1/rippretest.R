#PREPARE
R.Version()#for referencing, shows you which R version you are using
rm(list=ls())#removes any other items in your workspace
ls()#check whether workspace is empty


pretest<-read.csv("ripchild_pretest.csv", header=T)
#all the libraries necessary to run the analyses
library(lme4)
library(readr)
library(tidyverse)
library(sjPlot)
#install.packages('TMB', type='source')
library(ggthemes)
library(gridExtra)
library(reshape2)
library(car)
#library(psych)
library("ggpubr")
library(dplyr)


#Overview

str(pretest)

#recode variables
pretest$id <- as.factor(pretest$id)
pretest$sex <- as.factor(pretest$sex)
pretest$order <- as.factor(pretest$order)
str(pretest)

#aggregate data

pretest_ind <- pretest %>%
  group_by(id, sex, agemonths, order) %>%
  summarize(correct = mean(correct))

#Sex distribution

table(pretest_ind$sex)
table(pretest_ind$sex, pretest_ind$order)
boys <- subset(pretest_ind, sex == "m")
girls <- subset(pretest_ind, sex == "f")
boys_sticker <- subset(boys, order == "stickerfirst")
boys_pen <- subset(boys, order == "penfirst")
girls_sticker <- subset(girls, order == "stickerfirst")
girls_pen <- subset(girls, order == "penfirst")


#table S1 redone
summary(boys_sticker)
sd(boys_sticker$agemonths)
summary(boys_pen)
sd(boys_pen$agemonths)

summary(girls_sticker)
sd(girls_sticker$agemonths)
summary(girls_pen)
sd(girls_pen$agemonths)