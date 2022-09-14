# Preprocessing of annual surveys for CRAMPP
# Matt Kmiecik
# Started 14 SEPTEMBER 2022

source("r-prep.R") # Prepares R workspace

# Loads data - - - -

# parsed API data
load("../output/arm1_annuals_parsed.RData") # arm 1 annuals
load("../output/arm1_year1_followup_parsed.RData") # arm 1 non-standard y1 annuals
load("../output/arm2_annuals_parsed.RData") # arm 2 annuals
load("../output/short_annuals_parsed.RData") # shortened annuals

# crampp participant codes
crampp_codes <- read_excel("../data/crampp-codes.xlsx", sheet = "usethis")

# variables of interest across assessment visits and annual questionnaires were
# summarized in an excel document here:
data_dict <- read_excel(path = "../doc/annuals-compressed-data-dictionary.xlsx")

# converting into long format for easy data extraction
data_dict_long <- 
  data_dict %>% 
  select(questionid:shortened_annual_var) %>%
  pivot_longer(-questionid, names_to = "proj", values_to = "question") %>%
  left_join(
    ., 
    select(data_dict, questionid, custom_name, cat_1), 
    by = "questionid"
  ) %>%
  separate(cat_1, into = c("main", "sub"), remove = FALSE)

# Preprocessing - - - -

####################
#                  #
# Shortened Annual #
#                  #
####################

# To prevent duplicate values when doing pivot_wide later
this_dict <- data_dict_long %>% filter(proj == "shortened_annual_var")

# initial data processing
short_annual_data_wide <- 
  short_annuals_parsed %>% 
  filter(field_name %in% this_dict$question) %>% # narrows down
  mutate(year = as.numeric(regmatches(event_id, regexpr("\\d", event_id)))) %>% # retrieves year EVENT ID DOES NOT EXIST REPLACE WITH REDCAP_EVENT_NAME
  left_join(., icsiplus_data_dict_long, by = c("field_name" = "question")) %>%
  select(record, year, question = field_name, value, questionid:sub) %>%
  pivot_wider(id_cols = c(record, year), names_from = custom_name, values_from = value)







