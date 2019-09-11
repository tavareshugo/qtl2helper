
library(qtl2)

# Read example data and run scan
DOex <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip")
DOex_probs <- calc_genoprob(DOex)
DOex_scan1 <- scan1(DOex_probs, DOex$pheno)

# coefficients with and without SE
DOex_coefs <- scan1coef(DOex_probs[, 2], DOex$pheno[, 1, drop = FALSE])
DOex_coefs_se <- scan1coef(DOex_probs[, 2], DOex$pheno[, 1, drop = FALSE], se = TRUE)


test_that("coercion of scan1coef to tibble - no SE, no map",
          {

            coefs_tbl <- tidy(DOex_coefs)

            # expectations
            expect_true(tibble::is_tibble(coefs_tbl))
            expect_equal(nrow(coefs_tbl), nrow(DOex_coefs) * ncol(DOex_coefs))
            expect_equal(ncol(coefs_tbl), 3)
            expect_equal(colnames(coefs_tbl), c("marker", "coef", "estimate"))

            })

test_that("coercion of scan1coef to tibble - with SE, no map",
          {

            coefs_tbl <- tidy(DOex_coefs_se)

            # expectations
            expect_true(tibble::is_tibble(coefs_tbl))
            expect_equal(nrow(coefs_tbl), nrow(DOex_coefs) * ncol(DOex_coefs))
            expect_equal(ncol(coefs_tbl), 4)
            expect_equal(colnames(coefs_tbl), c("marker", "coef", "estimate", "SE"))

          })


test_that("coercion of scan1coef to tibble - no SE, with map",
          {

            coefs_tbl <- tidy(DOex_coefs, map = DOex$gmap)

            # expectations
            expect_true(tibble::is_tibble(coefs_tbl))
            expect_equal(nrow(coefs_tbl), nrow(DOex_coefs) * ncol(DOex_coefs))
            expect_equal(ncol(coefs_tbl), 5)
            expect_equal(colnames(coefs_tbl), c("marker", "chrom", "pos", "coef", "estimate"))

          })

test_that("coercion of scan1coef to tibble - with SE, with map",
          {

            coefs_tbl <- tidy(DOex_coefs_se, map = DOex$gmap)

            # expectations
            expect_true(tibble::is_tibble(coefs_tbl))
            expect_equal(nrow(coefs_tbl), nrow(DOex_coefs) * ncol(DOex_coefs))
            expect_equal(ncol(coefs_tbl), 6)
            expect_equal(colnames(coefs_tbl), c("marker", "chrom", "pos", "coef", "estimate", "SE"))

          })
