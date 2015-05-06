release_questions <- function() {
  c('Have you run    check(\'ltmle-dev\', cran = F)  ?', #doesn’t skip some tests
    'Have you run    run_examples(\'ltmle-dev\', test=F)   ?') #doesn’t skip some examples
}