# Ian Dworkin, October 24th 2011

# Here is a little script for Mike Payne, showing how I do lda using R. Here we are using the built in lda function in R (from the MASS library). If you would like the code for it in matrix algebra I can write that for you as well. 

# libraries
require(MASS) # contains lda() and qda() functions

# set working directory
setwd("~/Dropbox/WingBiometrics2014/")


# read in data. Here we are using 15 2d landmarks on Drosophila wings
wings <- read.csv("~/Dropbox/WingBiometrics2014/PredationWingsALL.csv", h=T)

# Some notes on this data.
# 1 - the landmarks were acquired "manually" from images using the ImageJ plug-in from Chris Klingenberg
# 2 - These images are from an analysis of phenotypic selection of juvenile mantids on Drosophila melanogaster.
# 3 - For the purposes of this overview, we only care about the shape data and sex of flies. We are going to ignore other variables (like survivorship)
# 4 - Side refers to Left, Right or unknown (from wings that fell off when the flies were being eaten by predator). We should take care of this to avoid pseudo-replication, but for the moment we will ignore it (but don't for a real analysis).
# 5 - Since we are starting by examining sex, we should also include size as variable, or at least use size corrected shape variables. Again, for the purposes of getting this all started, we are ignoring it.


# look at the structure of the data object.
str(wings)


# putting together right hand side of the model (to make life easier for inputting)
RHS <- paste("ProcCoord", 1:26, sep="", collapse= " + ")

# Putting the formula for the model as a model object (again to make life easier)
DiscrimEqn <- as.formula(paste("Sex ~ ", RHS, sep="", collapse=""))

# calling the function
linDiscrim.1 <- lda(DiscrimEqn, data = wings, method="mve")

#loadings for the linear discriminant.
linDiscrim.1$scaling


# set flag for cross validation
linDiscrim.2 <- lda(DiscrimEqn, data = wings, CV=T)

# calls for the classifications (not true sex of fly, but what lda is calling it)
linDiscrim.2$class

# Posterior probabilities for classifications
linDiscrim.2$posterior



# Which ones did it classify correctly
linDiscrim.2$class == wings$Sex

# Proportion correctly identified, ~ 0.936
sum(linDiscrim.2$class == wings$Sex)/length(wings$Sex)  

## Below is just Ian checking things, can ignore.

# Computing the discriminant function coefficient vector, a
# For the two group problem
# First compute means for males and females
males <- wings[wings$Sex=="M",9:34]
females <- wings[wings$Sex=="F",9:34]


male.means <- colMeans(males)
female.means <- colMeans(females)

# pooled covariance matrix
pooled.cov <- cov(wings[, 9:34]) # This is S
S.inv <- solve(pooled.cov) # If you set the tolerance too low, all hell breaks loose, so for this example I am using a covariance matrix of full rank,

# Now compute the discriminant function coefficient vector, a
# a = S[inverse]

a = S.inv %*% (male.means - female.means)

# rescaling it to compare to the original vector
a.prime = a*(34.488708/31.918649)
linDiscrim.1$scaling