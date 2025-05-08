bf(dna_per_ul ~ is_control + (1 | sample_id),
   shape ~ is_control + (1 | sample_id))

library(gamlss)
fit_gamlss <- gamlss(
  dna_per_ul ~ is_control + random(sample_id),
  sigma.formula = ~ is_control + random(sample_id),
  family = GA,   # Gamma distribution
  data = mutate(merged_dna_quants, across(where(is.character), as.factor))
)

summary(fit_gamlss)

library(emmeans)
emmeans(fit_gamlss, ~is_control)
predictAll(fit_gamlss)


get_brms_priors <- function(){
  #Fit FREQ Model
  fit_gamlss <- gamlss(
    dna_per_ul ~ re(fixed=~is_control, random=~1|sample_id),
    sigma.formula = ~ re(fixed=~is_control, random=~1|sample_id),
    family = GA,   # Gamma distribution
    data = mutate(merged_dna_quants, across(where(is.character), as.factor))
  )

  # Fixed effects for the location parameter (mu)
  mu_fixed <- summary(fit_gamlss$mu.coefSmo[[1]])$tTable
  mu_intercept <- c(mu_fixed["(Intercept)", 'Value'],
                    mu_fixed["(Intercept)", 'Std.Error'])
  mu_control <- c(mu_fixed["is_controlTRUE", 'Value'],
                    mu_fixed["is_controlTRUE", 'Std.Error'])
  
  # Fixed effects for the dispersion parameter (sigma/shape)
  sigma_fixed <- summary(fit_gamlss$sigma.coefSmo[[1]])$tTable
  sigma_intercept <- c(sigma_fixed["(Intercept)", 'Value'],
                       sigma_fixed["(Intercept)", 'Std.Error'])
  sigma_control <- c(sigma_fixed["is_controlTRUE", 'Value'],
                     sigma_fixed["is_controlTRUE", 'Std.Error'])
  
  #Random Effects 
  ranef(fit_gamlss$mu.coefSmo[[1]]) %>%
    apply(2, sd)
  
  ranef(fit_gamlss$sigma.coefSmo[[1]])
  
  coef(getSmo(fit_gamlss))
  ranef(getSmo(fit_gamlss))
  fixef(getSmo(fit_gamlss))
  
}