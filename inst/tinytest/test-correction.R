library("enrichwith")

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))

mod <- glm(conc ~ lot*log(u), data = clotting, family = Gamma)
emod <- enrich(mod, with = "all")

coefs <- coef(emod, model = "mean")
disp <- coef(emod, model = "dispersion")
coefs_bc <- c(coefs, disp) - get_bias_function(emod)(coefs, disp)
attributes(coefs_bc) <- NULL

mod_bc <- glm(conc ~ lot*log(u), data = clotting, family = Gamma, method = "brglmFit", type = "correction")

## bias corrected estimates computed using enrichwith are the same as those when having type = 'correction' in brglmFit
expect_equal(unname(coefs_bc), unname(coef(mod_bc, model = "full")), check.attributes = FALSE, tolerance = 1e-06)

