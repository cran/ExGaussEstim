y<-rnorm(400, 6.5, 1.5)
y_na<-c(NA, y[-1])
y_char<-c('a', y[-1])
y_factor<-sample.int(3, 400, replace = T)
y_small<-rnorm(15, 6.5, 1.5)
df<-data.frame(y, y_factor)


testthat::test_that("errors", {
  testthat::expect_error(
    QMLEEstim(y, "NEM"),
    "this is not a valid argument, must be 'NEMD', 'FMIN', or 'NLMI'"
  )
  testthat::expect_error(
    QMLEEstim(y_na,"NEMD"),
    "'y' must not contain NAs"
  )
  testthat::expect_error(
    QMLEEstim(y_char,"NEMD"),
    "'y' must be numeric"
  )
  
  testthat::expect_error(
    QMLEEstim(y_factor,"NEMD"),
    "'y' must be numeric"
  )
  testthat::expect_error(
    QMLEEstim(y_small,"NEMD"),
    "'y' contains less than 20 values"
  )
  testthat::expect_error(
    QMLEEstim(df,"NEMD"),
    "'y' must be a vector"
  )
  testthat::expect_error(
    BayesianExgaussian(400,y_na),
    "'x' contains NAs"
  )
  testthat::expect_error(
    BayesianExgaussian(400,y_char),
    "'x' must be numeric"
  )
  testthat::expect_error(
    BayesianExgaussian(400,y_factor),
    "'x' has less than 4 distinct values"
  )
  testthat::expect_error(
    BayesianExgaussian(15, y_small),
    "'x' contains less than 19 values"
  )
  testthat::expect_error(
    BayesianExgaussian(400,df),
    "'x' must be a vector"
  )
})


