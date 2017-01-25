#################################################
## Author: Hannah Gahagan and Nikhil Milind    ##
## Date: January 12 2016                       ##
#################################################

source("src/configuration.R")

set.seed(2014)
data("framing", package="mediation")

med.fit <- lm(emo ~ treat + age + educ + gender + income, data=framing)
out.fit <- glm(cong_mesg ~ emo + treat + age + educ + gender + income, data=framing, family=binomial("probit"))

med.out <- mediate(med.fit, out.fit, treat="treat", mediator="emo", robustSE=TRUE, sims=100) 
summary(med.out)

INS.10wk <- log(clinical(clin.name="INS.10wk", data.set=phenotypes))
GLU.4wk <- log(clinical(clin.name="GLU.4wk", data.set=phenotypes))
Fbn1 <- gene.exp(gene.name="Fbn1", data.set=adipose.rz)

model.data <- data.frame(INS.10wk, GLU.4wk, Fbn1)

med.fit <- lm(GLU.4wk ~ Fbn1, data=model.data)
out.fit <- glm(INS.10wk ~ Fbn1 + GLU.4wk, data=model.data)

med.out <- mediate(med.fit, out.fit, treat="Fbn1", mediator="GLU.4wk", robustSE=TRUE, sims=100)

summary(med.out)
