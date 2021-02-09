usethis::use_build_ignore("_devhistory.R")
usethis::use_git(message = ":tada: Initial commit")

usethis::use_git_ignore(".DS_Store")
usethis::use_build_ignore(".DS_Store")
usethis::use_git(message = ":see_no_evil: Ban .DS_Store files")

usethis::use_gpl3_license(name="Karine Princé")
usethis::use_git(message = ":page_facing_up: Add package license")

# documentation
usethis::edit_file("DESCRIPTION")
usethis::use_git(message = ":bulb: Edit package metadata")
usethis::use_package_doc()
devtools::document()
usethis::use_git(message = ":bulb: Update documentation")

usethis::use_readme_rmd()
usethis::use_vignette("tutorial")

# build package
devtools::check()
devtools::install()
devtools::build_readme()
usethis::use_r("fun_abundance_estimation.R")
usethis::use_r("fun_abundance_habitat_estimation.R")
usethis::use_git(message = "fun_abundance: start documentation")
usethis::use_test(name = "test_fun_abundance_estim")

devtools::load_all()

# CI
usethis::use_travis()
usethis::use_github_action_check_standard()
usethis::use_github_action("pkgdown")

# use pkgdown
pkgdown::build_site()
usethis::use_git_ignore("docs")
usethis::use_git_ignore("inst/docs")


