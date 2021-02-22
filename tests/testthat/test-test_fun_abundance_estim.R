context("Checking input data")


test_that("input data match expected format"{
  df%>%
    expect_s3_class("data.frame") %>%
    expect_length(5)

  expect_named(df, c("species", "site", "scheme","effort","N_obs"), ignore.order = TRUE, ignore.case = TRUE)
})
