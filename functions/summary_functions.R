create_map_age = function(age_max){
  # create map by 5-year age bands
  df_age_continuous <<- data.table(age_from = 0:age_max,
                                   age_to = 0:age_max,
                                   age_index = 0:age_max,
                                   age = 0:age_max)
  
  # create map for reporting age groups before 2020-09-02
  df_age_reporting_1 <<- data.table(age_from = c(0,5,15,25,35,45,55,65,75,85),
                                  age_to = c(4,14,24,34,44,54,64,74,84,age_max),
                                  age_index = 1:10,
                                  age_cat = c('0-4', '5-14', '15-24', '25-34', '35-44', '45-54', '55-64', '65-74', '75-84', '85+'))
  df_age_reporting_1[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age_cat"]
  df_age_reporting_1[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age_cat"]
  
  # create map for reporting age groups after 2020-09-02
  df_age_reporting_2 <<- data.table(age_from = c(0,5,10,18,30,50,65,75,85),
                                    age_to = c(4,9,17,29,49,64,74,84,age_max),
                                    age_index = 1:9,
                                    age_cat = c('0-4', '5-9', '10-17', '18-29', '30-49', '50-64', '65-74', '75-84', '85+'))
  df_age_reporting_2[, age_from_index := which(df_age_continuous$age_from == age_from), by = "age_cat"]
  df_age_reporting_2[, age_to_index := which(df_age_continuous$age_to == age_to), by = "age_cat"]
}
