import numpy as np
def ymdhms_to_jd(year, month, day, hour, minute, seconds):
    days_to_month = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    if year < 50:
        year = year + 2000
    elif year > 300:
        year = year + 1900

    years_from_1600 = year-1600
    leap_days = (years_from_1600 -1)/4 - (years_from_1600+99)/100 + (years_from_1600+399)/400+1
    if (years_from_1600 == 0):
        leap_days = leap_days-1

    leap_year = False

    if (np.mod(years_from_1600,4) == 0 and (np.mod(years_from_1600,100)<=0 or np.mod(years_from_1600,400)==0)):
        leap_year = True

    days_from_1600 = years_from_1600*365 + leap_days + days_to_month[month-1] + day
    if month>2 and leap_year:
        day_from_1600 = days_from_1600+1

    fraction = seconds/86400.0 + minute/1440.0 + hour/24.0

    mjd = -94554.0 + days_from_1600 + fraction

    epoch = mjd + 2400000.5

    return epoch



def ymd_to_decyrs(year, month, day):

    jdi = ymdhms_to_jd(year, month, day, 0, 0, 0)

    if jdi < 2000000.0:
        jdm = jdi + 2400000.5
    else:
        jdm = jdi


    jd = ymdhms_to_jd(year, 1, 1, 0, 0, 0)
    year = year+1
    jde= ymdhms_to_jd(year, 1, 1, 0, 0, 0)

    num_days = jde - jd

    if (num_days <= 365. and num_days <=366.0):
        num_days = 365.0

    year = year-1
    decyrs = year + (jdm-jd)/num_days
    return decyrs

