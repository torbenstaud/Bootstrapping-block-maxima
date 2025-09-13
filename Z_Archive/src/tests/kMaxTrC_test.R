# Load the testthat package
library(testthat)

# Test 1: Basic functionality with a regular input
test_that("kMaxTrC handles regular input correctly ;r = k", {
  sample <- c(1, 3, 2, 5, 6, 7, 8, 4)
  r <- 2
  k <- 2
  result <- kMaxTrC(sample, r, k)
  
  expected_result <- c(3, 3, 5, 5, 7, 8, 8, 6)
  expect_equal(result, expected_result)
})
# Test 2: Basic functionality with a regular input: r!= k
test_that("kMaxTrC handles regular input correctly ;r != k", {
  sample <- c(1, 3, 2, 5, 6, 7, 8, 4, 8, 12, 3, 1)
  r <- 3
  k <- 2
  result <- kMaxTrC(sample, r, k)
  
  expected_result <- c(3, 5, 6, 7, 7, 7, 8, 12, 12, 12, 8, 8)
  expect_equal(result, expected_result)
})



# Test 3: Edge case where n < k * r
test_that("kMaxTrC handles small input where n < k * r", {
  sample <- c(2, 9, 4, 8, 10)
  r <- 3
  k <- 2
  result <- kMaxTrC(sample, r, k)
  
  expected_result <- c(9, 9, 10, 10, 10)
  expect_equal(result, expected_result)
})

# Test 4: Case where the last block is smaller than r
test_that("kMaxTrC handles last block smaller than r", {
  sample <- c(5, 2, 8, 6, 4, 2, 4, 3)
  r <- 3
  k <- 2
  result <- kMaxTrC(sample, r, k)
  
  expected_result <- c(8, 8, 8, 6, 5, 5, 4, 4)
  expect_equal(result, expected_result)
})



# Run the tests
if(F){
  test_file("kMaxTrC_test.R")
}
