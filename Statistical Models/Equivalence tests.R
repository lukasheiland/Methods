##########################################################################################
# Equivvalence test ---------------------------------------------------------------------
##########################################################################################


# two packages
# library('TOSTER') # syntax
library('equivalence')



# Equivalence test ---------------------------------------------------------------------
# composite null hypothesis (= infinite null hypotheses)
# TOST := two one sided t-tests
data(ufc)

ht <- tost(ufc$Height.m.p, ufc$Height)
xyplot(ht)

ufc.ht <- ufc[!is.na(ufc$Height),]
equivalence.xyplot(ufc.ht$Height.m ~ ufc.ht$Height.m.p,
                   alpha=0.05, b0.ii=0.1, b1.ii=0.2,
                   xlab="Predicted height (m)",
                   ylab="Measured height (m)")
plot(ufc)



