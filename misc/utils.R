ensure_increasing_cumulative_deaths_origin = function(tmp)
{
  
  #
  # check if cumulative death is strictly decreasing from last date to first date
  for(Code in unique(tmp$State)){
    for(age_group in unique(tmp$Age.group)){
      
      tmp1 = subset(tmp, Age.group == age_group & State == Code)
      dates = unique(tmp1$date)
      
      for(t in 1:length(dates)){
        
        Date = rev(dates)[t]
        #print(Date)
        
        if(Date < max(dates)){
          
          if(is.na(tmp[Age.group == age_group & date == Date & State == Code, COVID.19.Deaths]))
            next
          if(is.na(tmp[Age.group == age_group & date == rev(dates)[t-1] & State == Code, COVID.19.Deaths]))
            next
          
          # if cumulative date at date t < date t - 1, fix cum deaths at date t - 1 to the one at date t.
          if(tmp[Age.group == age_group & date == Date & State == Code, COVID.19.Deaths] > tmp[Age.group == age_group & date == rev(dates)[t-1] & State == Code, COVID.19.Deaths]){
            tmp[Age.group == age_group & date == Date & State == Code,]$COVID.19.Deaths = tmp[Age.group == age_group & date == rev(dates)[t-1] & State == Code, COVID.19.Deaths]
          }
          
        }
        
      }
    }
  }
  
  
  return(tmp)
}

bugfix_nonincreasing_cumulative_deaths = function(tmp)
{
  # notices non stricly increasing cum deaths
  
  if(unique(tmp$Sex) == 'Female'){
    tmp[State == 'Oregon' & Age.group == '40-49 years' & date == as.Date('2021-03-13')]$COVID.19.Deaths = 14 
  }
  
  return(tmp)
}

fix_inconsistent_NA_between0 = function(tmp)
{
  
  tmp1 = tmp[, list(idx_last_0 = (which(COVID.19.Deaths == 0))), by = c("State", 'Age.group')]
  tmp1 = tmp1[, list(idx_last_0 = max(idx_last_0)), by = c("State", 'Age.group')]
  
  tmp2 = unique(select(tmp, "State", 'Age.group'))
  tmp1 = merge(tmp1, tmp2, by = c("State", 'Age.group'), all.y = T)
  
  tmp = merge(tmp, tmp1, by = c("State", 'Age.group'))
  tmp[, date_idx := 1:length(date), by = c("State", 'Age.group')]
  tmp[date_idx <= idx_last_0 & is.na(COVID.19.Deaths), COVID.19.Deaths := 0]
  
  tmp = select(tmp, -date_idx, -idx_last_0)
  
  return(tmp)
}

fix_inconsistent_NA_betweenpos = function(tmp)
{
  
  tmp1 = tmp[, list(idx_NA = (which(is.na(COVID.19.Deaths)))), by = c("State", 'Age.group')]
  tmp1 = tmp1[, list(idx_last_NA = max(idx_NA), idx_first_NA = min(idx_NA)), by = c("State", 'Age.group')]
  
  tmp2 = unique(select(tmp, "State", 'Age.group'))
  tmp1 = merge(tmp1, tmp2, by = c("State", 'Age.group'), all.y = T)
  
  tmp = merge(tmp, tmp1, by = c("State", 'Age.group'))
  tmp[, date_idx := 1:length(date), by = c("State", 'Age.group')]
  tmp[date_idx >= idx_first_NA & date_idx <= idx_last_NA & COVID.19.Deaths > 0 & COVID.19.Deaths < 15, COVID.19.Deaths := NA]
  
  stopifnot(unique(tmp[date_idx == (idx_first_NA - 1)]$COVID.19.Deaths)==0)
  
  tmp = select(tmp, -date_idx, -idx_first_NA, -idx_last_NA)
  
  return(tmp)
}

ensure_increasing_cumulative_deaths = function(tmp)
{
  
  dates = unique(tmp$date)
  
  for(t in 1:length(dates)){
    
    Date = rev(dates)[t]
    print(Date)
    #
    # check if cumulative death is strictly decreasing from last date to first date
    
    for(age_group in unique(tmp$age)){
      
      for(Code in unique(tmp$code)){
        
        if(Date < max(dates)){
          
          if(is.na(tmp[age == age_group & date == Date & code == Code, COVID.19.Deaths]))
            next
          if(is.na(tmp[age == age_group & date == rev(dates)[t-1] & code == Code, COVID.19.Deaths]))
            next
          
          # if cumulative date at date t < date t - 1, fix cum deaths at date t - 1 to the one at date t.
          if(tmp[age == age_group & date == Date & code == Code, COVID.19.Deaths] > tmp[age == age_group & date == rev(dates)[t-1] & code == Code, COVID.19.Deaths]){
            tmp[age == age_group & date == Date & code == Code,]$COVID.19.Deaths = tmp[age == age_group & date == rev(dates)[t-1] & code == Code, COVID.19.Deaths]
          }
          
        }
        
      }
    }
  }
  
  
  return(tmp)
}

map_statename_code = data.table(State = c(
  "Alabama"         ,         
  "Alaska"          ,         
  "Arizona"	        ,         
  "Arkansas"        ,         
  "California"      ,         
  "Colorado"        ,         
  "Connecticut"     ,         
  "Delaware"        ,         
  "Florida"         ,         
  "Georgia"         ,         
  "Hawaii"          ,         
  "Idaho"           ,         
  "Illinois"        ,         
  "Indiana"         ,         
  "Iowa"            ,         
  "Kansas"          ,         
  "Kentucky"        ,         
  "Louisiana"	      ,         
  "Maine"           ,         
  "Maryland"        ,         
  "Massachusetts"   ,         
  "Michigan"        ,         
  "Minnesota"       ,         
  "Mississippi"     ,         
  "Missouri"        ,         
  "Montana"         ,         
  "Nebraska"        ,         
  "Nevada"          ,         
  "New Hampshire"   ,         
  "New Jersey"      ,         
  "New Mexico"      ,         
  "New York"        ,
  "New York City"   ,
  "North Carolina"  ,         
  "North Dakota"    ,         
  "Ohio"            ,         
  "Oklahoma"        ,         
  "Oregon"	        ,         
  "Pennsylvania"    ,         
  "Rhode Island"    ,         
  "South Carolina"  ,         
  "South Dakota"    ,         
  "Tennessee"	      ,         
  "Texas"	          ,         
  "Utah"	          ,         
  "Vermont"         ,         
  "Virginia"        ,         
  "Washington"      ,         
  "West Virginia"   ,         
  "Wisconsin",	               
  "Wyoming"), 
  code = c(
    "AL",
    "AK",
    "AZ",
    "AR",
    "CA",
    "CO",
    "CT",
    "DE",
    "FL",
    "GA",
    "HI",
    "ID",
    "IL",
    "IN",
    "IA",
    "KS",
    "KY",
    "LA",
    "ME",
    "MD",
    "MA",
    "MI",
    "MN",
    "MS",
    "MO",
    "MT",
    "NE",
    "NV",
    "NH",
    "NJ",
    "NM",
    "NY",
    "NYC",
    "NC",
    "ND",
    "OH",
    "OK",
    "OR",
    "PA",
    "RI",
    "SC",
    "SD",
    "TN",
    "TX",
    "UT",
    "VT",
    "VA",
    "WA",
    "WV",
    "WI",
    "WY"))
