context("paretoconv")

test_that ("sanity checks", {
               expect_error (paretoconv (), "x must be supplied")
               expect_error (paretoconv (x=1), "a must be supplied")
               expect_error (paretoconv (x=1, a=1), "n must be supplied")

               # These errors are ordered (x, a, n)
               expect_error (paretoconv (n=1), "x must be supplied")

               expect_error (paretoconv (x="1", a=1, n=1), "x must be numeric")
               expect_error (paretoconv (x=1, a="1", n=1), "a must be numeric")
               expect_error (paretoconv (x=1, a=1, n="1"), "n must be numeric")
               expect_error (paretoconv (x=1, a=1:2, n=1), 
                             "a must be a single value only")
               expect_error (paretoconv (x=1, a=1, n=1:2), 
                             "n must be a single value only")
               expect_error (paretoconv (x=1, a=1:2, n=1:2), 
                             "a must be a single value only") # also ordered

               expect_error (paretoconv (x=1, a=1, n=1.5), 
                             "n must be an integer")
               expect_error (paretoconv (x=1, a=1, n=0), 
                             "n must be an integer > 0")
               expect_error (paretoconv (x=-1, a=1, n=1), "x must be positive")
               expect_error (paretoconv (x=-1, a=1, n=-1), "x must be positive")
})

test_that("ramsay-int", {
              expect_that (paretoconv (x=1, a=1, n=1), is.numeric)
              expect_that (paretoconv (x=1, a=0.5, n=1), is.numeric)
              expect_that (paretoconv (x=1, a=1, n=1, cdf=TRUE), is.numeric)
              expect_that (paretoconv (x=1, a=0.5, n=1, cdf=TRUE), is.numeric)
              expect_true (length (paretoconv (x=1:3, a=1, n=1)) == 3)
})
