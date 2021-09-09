context("test print methods")

## source(system.file("inst", "brglm0/brglm0.R", package = "brglm2"))
data("lizards", package = "brglm2")


types <- list("ML" = "(maximum likelihood)",
               "correction" = "(bias correction)",
               "AS_mean" = "(mean bias-reducing adjusted score equations)",
               "AS_median" = "(median bias-reducing adjusted score equations)",
               "AS_mixed" = "(mixed bias-reducing adjusted score equations)",
               "MPL_Jeffreys" = "(maximum penalized likelihood with Jeffreys'-prior penalty)")


for (type in names(types)) {
    liz <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time, family = binomial(),
               data = lizards,
               method = "brglmFit", type = type)


}
