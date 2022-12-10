
#turn appropriate columns into numeric that shouldn't be labeled as character
d$unit_number <- as.numeric(d$unit_number)
d$unit_number[is.na(d$unit_number)] <- -99
d$salinity <- as.numeric(d$salinity)
d$pHat25 <- as.numeric(d$pHat25)


#calc DIC in umol/kg
d$dic <- carb(flag = 8, d$pHat25, d$alk/1000000, T = d$spec_temp, S = d$salinity)$DIC * 1000000
#get logical vector of rows with data to calc in situ pH
isOK <- !is.na(d$alk) & !is.na(d$dic) & !is.na(d$insituTemp) & !is.na(d$salinity)
#calc insitu pH
if(any(isOK)){
  d$pHinsitu[isOK] <- carb(flag = 15, d$alk[isOK]/1000000, d$dic[isOK]/1000000, T = d$insituTemp[isOK], S = d$salinity[isOK])$pH
}
values$specData <- d