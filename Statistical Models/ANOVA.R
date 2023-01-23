## One-way ANOVA

?anova()

n <- 100
ngroups <- 3
N <- ngroups * n

x1 <- rnorm(n, 1, 1)
x2 <- rnorm(n, 2, 1)
x3 <- rnorm(n, 0, 1)

X <- data.frame(X = c(x1, x2, x3),
                Group = rep(c("x1", "x2", "x3"), each = n ))

mean1 <- mean(x1)
mean2 <- mean(x2)
mean3 <- mean(x3)

meangesamt <- mean(X$X)

# Gesamtvarianz
# var = sum(means - X)-
sb1 <- (mean1 - meangesamt)^2 * n +
       (mean2 - meangesamt)^2 * n +
       (mean3 - meangesamt)^2 * n
df1 <- (3-1) 
sb1 <- sb1/df # df
# durchschnittliche gruppenvarianz
df2 <- 3*(100-1)
sb2 <- sum((x1 - mean1)^2 + (x2 - mean2)^2 + (x3 - mean3)^2)/df2

# Gesamtvarianz durch durchschnittliche Gruppenvarianz
r <- df(sb1/sb2, df1, df2)
qf(0.95, df1, df2)

p <- pf(r, 2, 297)

# two way ANOVA?
