# API Calls for CRAMPP
# Matt Kmiecik
# Started: 13 SEPTEMBER 2022

# Purpose: run this script to refresh Redcap annual data

# Loads packages - - - -
source("r-prep.R") # Prepares R workspace

# API setup - - - -
api_url <- "https://survey.northshore.org/api/" # Redcap URL

# Tokens
arm1_api_token <- 
  read_lines(
    file = "../private/EH13-094-OCP-GRANT-api-token.txt" # token
  )

arm2_api_token <- 
  read_lines(
    file = "../private/EH13-094-Arm-2-api-token.txt" # token
  )

shortened_annual_api_token <- 
  read_lines(
    file = "../private/EH13-094-Shortened-Annual-api-token.txt" # token
  )

# Function for API calls - - - -
# PROBLEM IS HERE: https://education.arcus.chop.edu/redcap-r-windows/

#########################################
#                                       #
# Grabs annuals from EH13-094-OCP-GRANT #
#                                       #
#########################################

arm1_annuals <- 
  POST(
    url = api_url,
    body = list(
      token = arm1_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      "forms[0]" = "annual_questionnaire_combined060416", # form
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
arm1_annuals_parsed <- content(arm1_annuals, as = "parsed")

# Trying to get time stamps of assessment visit
# pulls all data bc why not
arm1_avisit1_api <- 
  POST(
    url = api_url,
    body = list(
      token = arm1_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      #"forms[0]" = "ic_problem_and_symptoms_indices",
      "events[0]" = "assessment_visit_1_arm_1",
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
arm1_avisit1_api_parsed <- content(arm1_avisit1_api, as = "parsed")

# SCREEN VISIT DATA
arm1_screenvisit_api <- 
  POST(
    url = api_url,
    body = list(
      token = arm1_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      "events[0]" = "screen_visit_arm_1",
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
arm1_screenvisit_api_parsed <- content(arm1_screenvisit_api, as = "parsed")

# NON-STANDARD FOLLOWUP
# Some participants were given original questionnaires at year 1 followup
# these variables do not start with an "a"
arm1_year1_followup <- 
  POST(
    url = api_url,
    body = list(
      token = arm1_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      "events[0]" = "year_1_followup_do_arm_1",  
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
arm1_year1_followup_parsed <- content(arm1_year1_followup, as = "parsed")


#########################################
#                                       #
# Grabs annuals from EH13-094- Arm 2    #
#                                       #
#########################################

arm2_annuals <- 
  POST(
    url = api_url,
    body = list(
      token = arm2_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      "forms[0]" = "annual_questionnaire_combined060416", # form annual_questionnaire_combined060416
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
arm2_annuals_parsed <- content(arm2_annuals, as = "parsed")

# Trying to get time stamps of assessment visit
# pulls all avisit 1 data bc why not
arm2_avisit1_api <- 
  POST(
    url = api_url,
    body = list(
      token = arm2_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      #"forms[0]" = "ic_problem_and_symptoms_indices",
      "events[0]" = "assessment_visit_1_arm_1",
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
arm2_avisit1_api_parsed <- content(arm2_avisit1_api, as = "parsed")

# SCREEN VISIT DATA
arm2_screenvisit_api <- 
  POST(
    url = api_url,
    body = list(
      token = arm2_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      #"forms[0]" = "ic_problem_and_symptoms_indices",
      "events[0]" = "screen_visit_arm_1",
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
arm2_screenvisit_api_parsed <- content(arm2_screenvisit_api, as = "parsed")


#######################################
#                                     #
# EH13-094 Shortened Annual Follow Up #
#                                     #
#######################################

short_annuals <- 
  POST(
    url = api_url,
    body = list(
      token = shortened_annual_api_token,
      content = "record",
      format = "csv",
      type = "eav", # for longitudinal output
      "forms[0]" = "shortened_annual_questionnaire", # 
      rawOrLabel = "raw",
      rawOrLabelHeaders = "raw",
      exportCheckboxLabel = "false",
      exportSurveyFields = "true", # reports timestamps etc.
      exportDataAccessGroups = "false",
      returnFormat = "csv"
    )
  )
# parses .$content
short_annuals_parsed <- content(short_annuals, as = "parsed")

# Saves out API data - - - -

# ARM 1 ANNUAL DATA
save(arm1_annuals_parsed, file = "../output/arm1_annuals_parsed.RData")
# ARM 1 ASSESSMENT VISIT 1 DATA 
save(arm1_avisit1_api_parsed, file = "../output/arm1_avisit1_api_parsed.RData")
# ARM 1 SCREEN VISIT DATA
save(arm1_screenvisit_api_parsed, file = "../output/arm1_screenvisit_api_parsed.RData")
# ARM 1 Year 1 FOLLOW UP AS ASSESSMENT VISIT 1
save(arm1_year1_followup_parsed, file = "../output/arm1_year1_followup_parsed.RData")

# ARM 2 ANNUAL DATA
save(arm2_annuals_parsed, file = "../output/arm2_annuals_parsed.RData")
# ARM 2 ASSESSMENT VISIT 1 DATA
save(arm2_avisit1_api_parsed, file = "../output/arm2_avisit1_api_parsed.RData")
# ARM 2 SCREEN VISIT DATA
save(arm2_screenvisit_api_parsed, file = "../output/arm2_screenvisit_api_parsed.RData")

# SHORTENED ANNUALS DATA
save(short_annuals_parsed, file = "../output/short_annuals_parsed.RData")

# Cleaning script objects - - - -
rm(
  arm1_annuals,
  arm1_annuals_parsed,
  arm1_avisit1_api,
  arm1_avisit1_api_parsed,
  arm2_annuals,
  arm2_annuals_parsed,
  arm2_avisit1_api,
  arm2_avisit1_api_parsed,
  short_annuals,
  short_annuals_parsed,
  api_url,
  arm1_api_token,
  arm2_api_token,
  shortened_annual_api_token,
  arm1_screenvisit_api,
  arm1_screenvisit_api_parsed,
  arm2_screenvisit_api,
  arm2_screenvisit_api_parsed,
  arm1_year1_followup,
  arm1_year1_followup_parsed
)