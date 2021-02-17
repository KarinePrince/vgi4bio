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
usethis::use_readme_rmd()
usethis::use_vignette("tutorial")
usethis::use_git(message = ":bulb: Update documentation")

rmarkdown::render("README.Rmd")     #Each time I edit the README.Rmd I have to update the .md with this command
usethis::use_build_ignore("README.html")
usethis::use_git_ignore("README.html")
usethis::use_git(message = ":pencil: Edit README")
system("git push")

# build package
devtools::check()
devtools::install()
devtools::build_readme()
usethis::use_r("fun_abundance_estimation.R")
usethis::use_git(message = "fun_abundance: start documentation")
usethis::use_git(message = "fun_abundance: update function")
usethis::use_test(name = "test_fun_abundance_estim")

usethis::use_r("fun_abundance_habitat_estimation.R")
usethis::use_git(message = "fun_abundance_hab: start documentation")
usethis::use_git(message = "fun_abundance_hab: update function")
usethis::use_test(name = "test_fun_abundance_hab_estim")

devtools::test()
# check package
devtools::load_all()

# CI
#Instead of locally (and each time we make changes) checking our  package with devtools::check() and running unit tests with devtools::test(), we can delegate these tasks to GitHub with the functions of the family: usethis::use_github_action_check_*()
usethis::use_travis()
usethis::use_github_action_check_standard()
usethis::use_github_action("pkgdown")

# use pkgdown
usethis::use_pkgdown()
pkgdown::build_site()
usethis::use_git_ignore("docs")
usethis::use_git_ignore("inst/docs")


