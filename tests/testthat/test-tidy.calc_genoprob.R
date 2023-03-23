
library(qtl2)

# Read example data and run scan
DOex <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/DOex.zip")
DOex_probs <- calc_genoprob(DOex)

# number of markers * samples * genotypes
n <- sum(unlist(lapply(DOex_probs, function(i) prod(dim(i)))))

# marker names
markers <- unlist(lapply(DOex_probs, function(i) dimnames(i)[[3]]))

test_that("coercion of calc_genoprob to tibble - no map",
          {

            probs_tbl <- tidy(DOex_probs)

            # expectations
            expect_true(tibble::is_tibble(probs_tbl))
            expect_true(all(probs_tbl$marker %in% markers))
            expect_equal(nrow(probs_tbl), n)
            expect_equal(ncol(probs_tbl), 4)
            expect_equal(colnames(probs_tbl), c("marker", "id", "genotype", "probability"))

            })

test_that("coercion of calc_genoprob to tibble - with map",
          {

            probs_tbl <- tidy(DOex_probs, map = DOex$gmap)

            # expectations
            expect_true(tibble::is_tibble(probs_tbl))
            expect_equal(nrow(probs_tbl), n)
            expect_equal(ncol(probs_tbl), 6)
            expect_equal(colnames(probs_tbl), c("marker", "chrom", "pos", "id", "genotype", "probability"))

          })
