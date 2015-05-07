release_questions <- function() {
  c('Have you run    check(\'ltmle-dev\', cran = F)  ?', #do not skip any tests
    'Have you run    run_examples(\'ltmle-dev\', test=F)   ?', #do not skip any examples
    'Have you checked for nonASCII characters? - see R packages Evernote') 
}