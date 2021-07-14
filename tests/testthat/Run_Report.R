# file_coverage function
require(covr)
require(testthat)
require(devtools)
require(dplyr)

covr::file_coverage("R/calib_area.R",
                    test_files = "tests/testthat/test-calib_area.R") %>% report()

# To view if a misspell option as c("maks") send a warning and returns a SpatVector
covr::file_coverage("R/calib_area.R",
                    test_files = "tests/testthat/test0-calib_area.R") %>% report()

# To view if option as c("mcp") returns a SpatVector
covr::file_coverage("R/calib_area.R",
                    test_files = "tests/testthat/test1-calib_area.R") %>% report()

# To view if option "bmcp" returns a SpatialPolygon
covr::file_coverage("R/calib_area.R",
                    test_files = "tests/testthat/test2-calib_area.R") %>% report()

# To see if multiple options as c("mask", "clusters") returns an error
covr::file_coverage("R/calib_area.R",
                    test_files = "tests/testthat/test3-calib_area.R") %>% report()

# To see if c("mask", clusters, "clusters") returns a "Spatvector"
covr::file_coverage("R/calib_area.R",
                    test_files = "tests/testthat/test4-calib_area.R") %>% report()

# To see if we use a rgdal::readOGR object returns to "Spatvector"
covr::file_coverage("R/calib_area.R",
                    test_files = "tests/testthat/test5-calib_area.R") %>% report()



