


# potential arguments:
sdm_directory <- function(){
  fs::dir_create("./Project")
  fs::dir_create(here::here("Output/Intermediate"))
  fs::dir_create(here::here("Output/Report"))
  fs::dir_create(here::here("Output/Figures"))
}
