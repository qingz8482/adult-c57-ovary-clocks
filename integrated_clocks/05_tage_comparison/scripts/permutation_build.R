# Use the available system compiler rather than the unavailable conda compiler name.
compiler<-Sys.which('g++');stopifnot(nzchar(compiler))
makevars<-tempfile('compare_tage_makevars_',tmpdir='/tmp')
writeLines(paste(c('CXX','CXX11','CXX14','CXX17','CXX20','SHLIB_CXXLD'),compiler,sep=' = '),makevars)
Sys.setenv(R_MAKEVARS_USER=makevars)
