# Wrapper function to invoke "helloA" at the shell.
HelloWorld <- function() {
  system(paste(getwd(),"code/HelloWorld", sep="/"))
}