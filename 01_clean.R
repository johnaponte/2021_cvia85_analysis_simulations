#' ---
#' title: Clean directories
#' author: Created by repana::makestructure()
#' date:  2021-11-09 
#' ---
#' 
#' Clean the directories included in the
#' __clean_before_new_analysis__ section of config.yml
#' 
#' Execution date:
{{ format(Sys.time()) }}
#' 
repana::clean_structure()
