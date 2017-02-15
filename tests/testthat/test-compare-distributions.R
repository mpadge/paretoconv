context('compare-distributions')

test_that ('sanity checks', {
    expect_error (pdf_integral (), 'Value of a must be given')
    expect_error (compare_conv_distributions (), 
                  'Two distributions must be provided')
})

test_that ('compare_conv_distributions', {
    data ('native_american', package='poweRlaw')
    m1 <- poweRlaw::displ$new (native_american$Cas)
    m1$setXmin (3)
    m1$setPars (poweRlaw::estimate_pars (m1))
    x <- sort (unique (native_american$Cas))
    y <- paretoconv (x=x, a=m1$getPars (), x0=3, n=0, cdf=FALSE)
    #compare_conv_distributions (m1, y [1:(length (y) - 1)])

    expect_error (compare_conv_distributions (y, m1), 
                  'First argument must be a poweRlaw displ object')
    expect_error (compare_conv_distributions (m1, y [1:(length (y) - 1)]),
                  'Length of second argument does not correspond to data in first argument')
    expect_is (compare_conv_distributions (m1, y), 'compare_distributions')
})

test_that ('pdf-integrals', {
    res <- pdf_integral (a=2.5, x0=7, n=0, tol=1e-3)
    expect_is (res, 'numeric')
    expect_length (res, 1)
})
