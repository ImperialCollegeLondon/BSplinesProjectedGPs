create_map_age = function(age_max){
  # create map by 5-year age bands
  df_age_continuous <<- data.table(age_from = 0:age_max,
                                   age_to = 0:age_max,
                                   age_index = 0:age_max,
                                   age = 0:age_max)
  
  # create map for reporting age groups before 2020-09-02
  df_age_reporting <<- data.table(age_from = c(0,1, 5,15,25,35,45,55,65,75,85),
                                  age_to = c(0,4,14,24,34,44,54,64,74,84,age_max),
                                  age_index = 1:11,
                                  age_cat = c('0-0', '1-4', '5-14', '15-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '85+'))
  df_age_reporting[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age_cat"]
  df_age_reporting[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age_cat"]
  
}
