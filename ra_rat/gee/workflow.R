#!/usr/bin/Rscript
library(here)
library(stringr)

gee_script      <- here("gee", "gee_by_group.R")

mgs_profile     <- here("data", "mgs.profile")
ko_profile      <- here("data", "ko.profile")
module_profile  <- here("data", "module.profile")
pathway_profile <- here("data", "pathway.profile")

profile <- c(mgs_profile, ko_profile, module_profile, pathway_profile)
level <- c("mgs", "ko", "module", "pathway")

i = 0
for(level_ in level)
{
  i <- i + 1 
  gee_str <- str_c("Rscript ", gee_script,
                  " --profile ", profile[i],
                  " --level ", level_,
                  " --outdir ", here("gee", str_c(level_, "_marker_by_group/")))
  # print(gee_str)
  system(gee_str)
}


