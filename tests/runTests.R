## unit tests will not be done if testthat is not available
if(require("testthat", quietly=TRUE)) {
    test_dir("inst/unitTests/", reporter="minimal")
    library(hiAnnotator)    
    test_package("hiAnnotator")    
} else {
    warning("cannot run unit tests -- package 'testthat' is not available")
}
