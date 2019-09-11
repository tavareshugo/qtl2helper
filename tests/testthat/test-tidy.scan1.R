
library(qtl2)

# Read example data and run scan
DOex <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip")
DOex_probs <- calc_genoprob(DOex)
DOex_scan1 <- scan1(DOex_probs, DOex$pheno)


test_that("coercion of scan1 to tibble - no map",
          {

            scan1_tbl <- tidy(DOex_scan1)

            # expectations
            expect_true(tibble::is_tibble(scan1_tbl))
            expect_equal(nrow(scan1_tbl), nrow(DOex_scan1) * ncol(DOex_scan1))
            expect_equal(ncol(scan1_tbl), 3)
            expect_equal(colnames(scan1_tbl), c("marker", "pheno", "LOD"))

            })

test_that("coercion of scan1 to tibble - add gmap",
          {

            scan1_tbl <- tidy(DOex_scan1, map = DOex$gmap)

            # expectations
            expect_true(tibble::is_tibble(scan1_tbl))
            expect_equal(nrow(scan1_tbl), nrow(DOex_scan1) * ncol(DOex_scan1))
            expect_equal(ncol(scan1_tbl), 5)
            expect_equal(colnames(scan1_tbl), c("marker", "chrom", "pos", "pheno", "LOD"))

          })


