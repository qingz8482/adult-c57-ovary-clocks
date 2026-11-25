source(file.path(dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))),'common.R'))
source(file.path(here,'reference_setup.R'))
for(s in sort(list.files(file.path(here,'reference_renderers'),pattern='[.]R$',full.names=TRUE))){if(as.integer(substr(basename(s),1,2))>as.integer(Sys.getenv('REFERENCE_END', '99')))next;if(as.integer(substr(basename(s),1,2))<as.integer(Sys.getenv('REFERENCE_START', '1')))next;message('Rendering original reference layout: ',basename(s));source(s,local=new.env(parent=.GlobalEnv))}
message('PASS: reference-layout rendering complete.')
