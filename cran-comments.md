## Resubmission

This is a resubmission. Regarding the three mentioned issues:

* The main methodological reference for the package is not ready yet, so I have not added any Reference field (in DESCRIPTION) yet. My understanding is that it is not a mandatory field.

* I have replaced 'glm' by glm() in the Description field of DESCRIPTION file.

* Replaced the calls to `stats:::vcov.glm(resiso)` by `summary(resiso)$cov.scaled` (in man/confint.cirls.Rd) to avoid calling unexported objects.


## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a first submission
