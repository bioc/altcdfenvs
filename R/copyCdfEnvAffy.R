copyCdfEnvAffy <- function(acdfenv) {

  r <- acdfenv
  r@envir <- copyEnv(acdfenv@envir)

  return(r)
}
