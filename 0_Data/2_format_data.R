l_models <- list.files("/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/")
l_hist <- list.files("./", pattern = "historical") 
l_rcp85 <- list.files("./", pattern = "rcp85") 
hist_split <- strsplit(l_hist, split="_")
rcp85_split <- strsplit(l_rcp85, split="_")
l_run_names <- lapply(hist_split,'[',5)


find_file <- function(run_name, l_files){
  indice <- lapply(l_files,function(x)grep(x, pattern=run_name))
  indice <- which(unlist(lapply(indice,function(x) length(x) > 0))) 
  l_files[indice]
}

merge_files <- function(file1, file2, output){
  system(paste("cdo mergetime",file1, file2, output))
}

merge_run <- function(run_name, l_hist, l_rcp85){
  f_hist <- find_file(run_name, l_hist)
  f_rcp85 <- find_file(run_name, l_rcp85)
  stopifnot(length(f_hist) == 1 & length(f_rcp85) == 1)
  output <- paste(run_name, ".nc", sep="")
  merge_files(f_hist, f_rcp85, output)
  output
}

seldate <- function(input, output, start="1850-01-01T00:00:00", end="2100-12-31T23:59:59"){
  cmd <- (paste("cdo seldate,", start, ",",  end, " ",  input, " ", output, sep=""))
  print(cmd)
  system(cmd)
  output
}

format_run <- function(run_name, l_hist, l_rcp85){
  seldate(merge_run(run_name, l_hist, l_rcp85), paste(run_name, "_ok.nc", sep=""))
}

curf <- getwd()
l_models <- list.files("/mnt/nfs/d1/vdr/DATA/CMIP5/Origin/")
for (f in l_models){
  print(f)
  setwd(f)
  l_hist <- list.files("./", pattern = "historical") 
  l_rcp85 <- list.files("./", pattern = "rcp85") 
  hist_split <- strsplit(l_hist, split="_")
  rcp85_split <- strsplit(l_rcp85, split="_")
  l_run_names <- lapply(hist_split,'[',5)
  try(
      l_formated <- lapply(l_run_names, format_run, l_hist=l_hist, l_rcp85=l_rcp85)
  )
  setwd(curf)
}
