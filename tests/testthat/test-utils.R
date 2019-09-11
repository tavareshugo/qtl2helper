library(qtl2)

# Read example data and run scan
DOex <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip")
DOex_probs <- calc_genoprob(DOex)
DOex_scan1 <- scan1(DOex_probs, DOex$pheno)

scan1_tbl <- tidy(DOex_scan1)


test_that("join data.frame and list of genetic map positions with .join_map_by_marker",
          {
            scan1_tbl_map <- .join_map_by_marker(scan1_tbl, DOex$gmap)

            expect_true(tibble::is_tibble(scan1_tbl_map))
            expect_warning(.join_map_by_marker(scan1_tbl, DOex$gmap[2]))
            expect_error(.join_map_by_marker(scan1_tbl, DOex$gmap$`2`))

          })
