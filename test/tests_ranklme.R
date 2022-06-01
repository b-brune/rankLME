# tests for ranklme function
library(rankLME)

set.seed(99099)
d <- raw_data(400, 4, 2, 20, 20)

r <- ranklme(d$X, d$y, d$Z, d$g)

r2 <- ranklme(d$X, d$y, d$Z, d$g, weighted = TRUE)

set.seed(1213)
new_order <- sample(1:400)
d_perm <- lapply(d[c(1:3, 5)], function(x) {
  if (is.matrix(x)) x[new_order, , drop=FALSE] else x[new_order]
  })

r3 <- ranklme(d_perm$X, d_perm$y, d_perm$Z, d_perm$g) 

new_order <- c(2, 1, 3:400)
d_perm2 <- lapply(d[c(1:3, 5)], function(x) {
  if (is.matrix(x)) x[new_order, , drop=FALSE] else x[new_order]
})

r4 <- ranklme(d_perm2$X, d_perm2$y, d_perm2$Z, d_perm2$g) 

# Rankbased regression is permutation invariant which makes sense. But why isn't 
# the ranklme version of the function?
rfit(y ~ X[, -1], d)
rfit(y ~ X[, -1], d_perm)
ranklm(d$X, d$y, mean_function="median")$coef
ranklm(d_perm$X, d_perm$y, mean_function="median")$coef
