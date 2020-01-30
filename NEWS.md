# sticRs 1.1.3

* Bug fix when using set_param for associated plant

# sticRs 1.1.2

* optimi_stics now can restart automatically with new start values if the optimization process returns an error (e.g. singular matrix)  
* Add possibility to use integer parameter type in optimi_stics

# sticRs 1.1.1

* Fix several bugs in parameter optimization
* Best optimization is now returned by optimi_stics instead of the last one
* Add tests for optimization -> parameters min/max/start values

# sticRs 1.1.0
Add parameter optimization (adapted to intercropping)

# sticRs 1.0.2
Add vignettes for introduction and for sensitivity analyses.

# sticRs 1.0.1

Add exact match to set_* and read_param to avoid issue when using higher level functions such as stics_eval with similar parameter names (e.g. P_hauteur_threshold and P_hauteur_threshold_2).

# sticRs 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* First stable release of the package
