download_vaccine_data_state = function()
{
  data = read.csv(file = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/us_state_vaccinations.csv')
  write.csv(data, '~/git/covid19Vaccination/inst/data/us_state_vaccinations_210611.csv', row.names = F)
}