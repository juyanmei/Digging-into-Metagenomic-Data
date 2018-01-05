library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(stringr)
library(broom)
library(lubridate)
library(geepack)
library(openxlsx)
library(argparser)
library(here)

rm(list = ls())
p <-
  arg_parser("find marker using geeglm method for ra_rat project")
p <- add_argument(
  p,
  c("--profile", "--level", "--outdir"),
  help = c(
    "gene, mgs, ko, module, pathway profile",
    "gene, mgs, ko, module, pathway level",
    "output directory"
  )
)
argv <- parse_args(p)

profile <- argv$profile
level   <- argv$level
out_dir <- argv$outdir

# test
#level   <- "mgs"
#profile <- here("data", "mgs.profile")
#out_dir <- here("gee", "mgs_marker/")

if (!dir.exists(out_dir))
  dir.create(out_dir)

## loda data and tidy data
baseinfo_csv <- here("data", "rat_245S_baseinfo.csv")
footinfo_csv <- here("data", "rat_49S_footinfo.csv")

sample_info <- read_csv(baseinfo_csv)
foot_info   <-
  read_csv(footinfo_csv) %>%
  select(group, order_id, arthritis_score, day) %>%
  mutate(day = mdy(day)) %>%
  filter(
    day == ymd(20170720) |
      day == ymd(20170729) |
      day == ymd(20170804) |
      day == ymd(20170810) |
      day == ymd(20170819)
  ) %>%  # mdy("7/20/2017") == ymd(20170720) is TRUE
  mutate(seqtime_point =
           case_when(
             day == ymd(20170720) ~ 1,
             day == ymd(20170729) ~ 2,
             day == ymd(20170804) ~ 3,
             day == ymd(20170810) ~ 4,
             day == ymd(20170819) ~ 5
           )) %>%
  #mutate(arthritis_score = na_if(arthritis_score, 0))
  replace_na(list(arthritis_score = 0))

score_info <- right_join(sample_info, foot_info)
write_csv(score_info,
          here("data", "rat_245S_scoreinfo.csv"),
          col_names = TRUE)

# normalize score_info
score_info_norm <-
  mutate(score_info, arthritis_score = arthritis_score / sum(arthritis_score))

profile_info <-
  read_delim(profile, delim = "\t") %>%
  rename(marker_id = X1)

abun_norm <- function(x)
  x / sum(x)
profile_info_norm <-
  profile_info %>%
  mutate_at(vars(-marker_id), funs(mean))

# transform
profile_info_t <-
  profile_info %>%
  gather(fastq_id, marker_abun, -marker_id) %>%
  spread(marker_id, marker_abun)

profile_info_norm_t <-
  profile_info_norm %>%
  gather(fastq_id, marker_abun, -marker_id) %>%
  spread(marker_id, marker_abun)

score_profile_info <-
  left_join(score_info, profile_info_t) %>%
  mutate(group = factor(
    group,
    levels = c("Model", "HFD", "ZQFTN", "GJK", "MTX", "Lcasei", "Con")
  )) %>%
  arrange(group, order_id)
write_csv(
  score_profile_info,
  str_c(out_dir, "rat_245S_", level, "_score_profile_info.csv"),
  col_names = TRUE
)

score_profile_info_norm <-
  left_join(score_info_norm, profile_info_norm_t) %>%
  mutate(group =
           factor(
             group,
             levels = c("Model", "HFD", "ZQFTN", "GJK", "MTX", "Lcasei", "Con")
           )) %>%
  arrange(group, order_id)

no_zero <- function(x)
  # x[x == 0] <- 10 ^ (-20) # wrong
  x <- replace(x, x == 0, 10 ^ (-20))

score_profile_info_no_zero <-
  score_profile_info %>%
  mutate_at(
    vars(
      -sample_name,
      -sample_id,
      -library_id,
      -fastq_id,
      -order_id,
      -group,
      -seqtime_point,
      -arthritis_score,
      -day
    ),
    funs(no_zero)
  )

do_geeglm <- function(score_profile_, group_)
{
  count              <- 0
  all_gee_info_df    <- data_frame()
  marker_gee_info_df <- data_frame()
  marker_pvalue_df   <- data_frame()
  
  fit_arthritis <-
    formula(arthritis_score ~ marker_abun * group + day)
  
  for (i in 10:ncol(score_profile_))
  {
    groups <- c("Model", group_)
    score_abun <-
      score_profile_ %>%
      filter(group %in% groups) %>%
      select(seq(9), i) %>%
      mutate(group = factor(group, levels = groups)) # %>% View(.)
    #stop()
    marker_id <- colnames(score_abun)[10]
    colnames(score_abun)[10] <- "marker_abun"
    
    score_abun_summary <-
      score_abun %>%
      group_by(group) %>%
      summarise(count = sum(marker_abun > 0)) #%>% View(.)
    #stop()
    
    #if ((sum(score_abun_summary$count) > 10) &
    if ((sum(score_abun_summary$count) > length(score_abun$marker_abun) * 0.1) &
        (all(score_abun_summary$count > 0)))
    {
      fit_result <-
        tidy(
          geeglm(
            fit_arthritis,
            data = score_abun,
            id = order_id,
            family = gaussian,
            corstr = "ar1",
            scale.fix = TRUE
          )
        ) %>%
        mutate(term = ifelse(term == "marker_abun" , marker_id, term)) #%>% View(.)
      #stop()
      
      if ((filter(fit_result, term == marker_id)$p.value  < 0.05) &
          (filter(fit_result, term == str_c("group", group_))$p.value < 0.05) &
          (filter(fit_result, term == str_c("marker_abun:group", group_))$p.value < 0.05) &
          (log10(abs(
            filter(fit_result, term == marker_id)$estimate
          )) > 4))
      {
        marker_pvalue_df <-
          fit_result %>%
          filter(term == marker_id) %>%
          bind_rows(marker_pvalue_df, .)
        
        marker_gee_info_df <-
          fit_result %>%
          bind_rows(marker_gee_info_df, .)
      }
      
      all_gee_info_df <-
        fit_result %>%
        mutate(marker_id = marker_id) %>%
        bind_rows(all_gee_info_df, .)
      
      
      count <- count + 1
      print(str_c(count, " : ", i, " : ", marker_id))
    }
  }
  
  return(list(all_gee_info_df, marker_gee_info_df, marker_pvalue_df))
}

lca_gee <- do_geeglm(score_profile_info, "Lcasei")
con_gee <- do_geeglm(score_profile_info, "Con")
hfd_gee <- do_geeglm(score_profile_info, "HFD")
mtx_gee <- do_geeglm(score_profile_info, "MTX")
gjk_gee <- do_geeglm(score_profile_info, "GJK")
zqf_gee <- do_geeglm(score_profile_info, "ZQFTN")

write_csv(lca_gee[[1]], str_c(out_dir, "lca_all_gee_info.csv"), col_names = TRUE)
write_csv(con_gee[[1]], str_c(out_dir, "con_all_gee_info.csv"), col_names = TRUE)
write_csv(hfd_gee[[1]], str_c(out_dir, "hfd_all_gee_info.csv"), col_names = TRUE)
write_csv(mtx_gee[[1]], str_c(out_dir, "mtx_all_gee_info.csv"), col_names = TRUE)
write_csv(gjk_gee[[1]], str_c(out_dir, "gjk_all_gee_info.csv"), col_names = TRUE)
write_csv(zqf_gee[[1]], str_c(out_dir, "zqf_all_gee_info.csv"), col_names = TRUE)

write_csv(lca_gee[[2]],
          str_c(out_dir, "lca_marker_gee_info.csv"),
          col_names = TRUE)
write_csv(con_gee[[2]],
          str_c(out_dir, "con_marker_gee_info.csv"),
          col_names = TRUE)
write_csv(hfd_gee[[2]],
          str_c(out_dir, "hfd_marker_gee_info.csv"),
          col_names = TRUE)
write_csv(mtx_gee[[2]],
          str_c(out_dir, "mtx_marker_gee_info.csv"),
          col_names = TRUE)
write_csv(gjk_gee[[2]],
          str_c(out_dir, "gjk_marker_gee_info.csv"),
          col_names = TRUE)
write_csv(zqf_gee[[2]],
          str_c(out_dir, "zqf_marker_gee_info.csv"),
          col_names = TRUE)

write_csv(lca_gee[[3]],
          str_c(out_dir, "lca_marker_pvalue_info.csv"),
          col_names = TRUE)
write_csv(con_gee[[3]],
          str_c(out_dir, "con_marker_pvalue_info.csv"),
          col_names = TRUE)
write_csv(hfd_gee[[3]],
          str_c(out_dir, "hfd_marker_pvalue_info.csv"),
          col_names = TRUE)
write_csv(mtx_gee[[3]],
          str_c(out_dir, "mtx_marker_pvalue_info.csv"),
          col_names = TRUE)
write_csv(gjk_gee[[3]],
          str_c(out_dir, "gjk_marker_pvalue_info.csv"),
          col_names = TRUE)
write_csv(zqf_gee[[3]],
          str_c(out_dir, "zqf_marker_pvalue_info.csv"),
          col_names = TRUE)

write_lines(lca_gee[[3]]$term, str_c(out_dir, "Lcasei.txt"))
write_lines(con_gee[[3]]$term, str_c(out_dir, "Normal.txt"))
write_lines(hfd_gee[[3]]$term, str_c(out_dir, "HFD.txt"))
write_lines(mtx_gee[[3]]$term, str_c(out_dir, "MTX.txt"))
write_lines(gjk_gee[[3]]$term, str_c(out_dir, "GJK.txt"))
write_lines(zqf_gee[[3]]$term, str_c(out_dir, "ZQFTN.txt"))

save(list = ls(),
     file = str_c(out_dir, level, "_marker.Rdata"))
