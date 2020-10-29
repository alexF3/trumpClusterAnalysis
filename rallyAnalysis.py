# Alex Rich, MPH
# alex.rich@unc.edu
# updated 29 oct 20
# 
# Purpose: visualize spread of Covid-19, in this case near Trump rallies in WI 
# Inputs:  Johns Hopkins CSSE data (https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv',\
#               out = fpData +'covid_hopkins_overTime_CONFIRMED.csv)
#           US Census county shapefiles
#           US Census county populations and areas
#           Trump rally dates and addresses
#           API key for googlemaps geolocatoin
#           API key for openrouteservice drive time polygon API
#           (backup) API key for MapBox drive time polygon API

googleMapsKey = <enter gMaps key>
openrouteserviceKey = <enter openrouteservice key>

# Outputs:
# .html Leaflet map of counties
# .png of infection trends for counties that overlap a 1-hr drive polgygon vs. those with centroids within 150 miles
#   normalized to rolling-7-day-infection-avg for each category on the day of the rally in question
#   (e.g. if there avg-7-day-infection-count were 10 on the day of the rally and 20 by 14 days later, the end point would up "up by 100%")
# .csv of "local" counties with infection counts on key dates for Relative Risk calculations
# .csv of "surrounding" counties with infection counts on key dates for Relative Risk calculations
# .csv of summary findings including relative risk ratios before and after rally and relative population density (local / surrounding)

## Package Versions:
# wget==3.2
# pandas==1.1.1
# geopandas==0.8.1
# Shapely==1.7.1
# numpy==1.19.1
# geopy==1.21.0
# folium==0.9.1

import os
import wget
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import shape
from shapely.ops import cascaded_union
import numpy as np
import geopy.distance

import shapely
import matplotlib.pyplot as plt
from shapely.geometry import Point
import json
import datetime as dt
import folium

# ## Update Hopkins Cases and Deaths files (if needed)

fpData = os.getcwd()+'/data/'
# ## Get daily update of hopkins time series file for confirmed US cases
# wget.download('https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv',\
#               out = fpData +'covid_hopkins_overTime_CONFIRMED.csv')
# ## if file saves as a ____(01).csv, delete the old file and rename new to "covid_hopkins_overTime_CONFIRMED.csv"
# if os.path.exists(fpData + "covid_hopkins_overTime_CONFIRMED (1).csv"):
#     os.remove(fpData + "covid_hopkins_overTime_CONFIRMED.csv")
#     os.rename(fpData + "covid_hopkins_overTime_CONFIRMED (1).csv",fpData + "covid_hopkins_overTime_CONFIRMED.csv")
# fpData = os.getcwd()+'/data/'
# ## Get daily update of hopkins time series file for confirmed US cases
# wget.download('https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv',\
#               out = fpData +'covid_hopkins_overTime_DEATHS.csv')
# ## if file saves as a ____(01).csv, delete the old file and rename new to "covid_hopkins_overTime_CONFIRMED.csv"
# if os.path.exists(fpData + "covid_hopkins_overTime_DEATHS (1).csv"):
#     os.remove(fpData + "covid_hopkins_overTime_DEATHS.csv")
#     os.rename(fpData + "covid_hopkins_overTime_DEATHS (1).csv",fpData + "covid_hopkins_overTime_DEATHS.csv")
      

# Read in new US Confirmed Timeseries data
covidData = pd.read_csv(fpData +'covid_hopkins_overTime_CONFIRMED.csv',dtype={'FIPS':str})
covidData.FIPS = covidData.FIPS.str.split('.').str[0].str.zfill(5)


covidDeaths = pd.read_csv(fpData +'covid_hopkins_overTime_DEATHS.csv',dtype={'FIPS':str})
covidDeaths.FIPS = covidDeaths.FIPS.str.split('.').str[0].str.zfill(5)

# collect dates from timeseries file
dates = []
for i in covidData.columns:
    if '/' in i:
        dates.append(i)
        
fipsList = covidData[~covidData.FIPS.isnull()].FIPS.unique().tolist()


## Note: in rare cases the "new case" count goes negative!  Believe this to be a function of corrections to previously erronious data

#get county_df from shapefile
fp = os.getcwd() + '/data/US Counties Shapefiles'
county_df = gpd.read_file(fp,dtype={'STATEFP':str,'COUNTYFP':str})
# check data type so we can see that this is not a normal dataframe, but a GEOdataframe
county_df.NAME = county_df.NAME.str.lower()
county_df['fips'] = county_df.STATEFP + county_df.COUNTYFP

# Get trump rally data
locs = pd.read_csv(os.getcwd()+'/data/trumpRallies.csv')
locs['day'] = locs.date.str.split('[').str[0]
locs['day'] = pd.to_datetime(locs.day)
locs['rallyDate'] = pd.to_datetime(locs['day'])
locs['address'] = locs.location + ', ' + locs.city + ', ' + locs.state
locs['filename'] = locs['city'] + '_' + locs.state

# Read in county population and area data from US Census
countyPops = pd.read_csv(os.getcwd() + '/data/co-est2019-alldata.csv',encoding='Latin1',\
                        dtype={'STATE':str,'COUNTY':str})
countyPops['fips'] = countyPops.STATE + countyPops.COUNTY

######################### NOTE: REMOVED THIS FOR SIMPLICITY...CAN BE ADDED BACK IN TO EXLCUDE COUNTIES "LOCAL" TO OTHER RALLIES

# Read in frame of all of the counties that are local to Trump rallies (as of 25 Oct 20)
# allLocalCounties = pd.read_csv(os.getcwd()+'/data/all_local_counties_FIPS.csv',dtype={'localRallyCounties':str})
# allLocalCounties.localRallyCounties = allLocalCounties.localRallyCounties.str.zfill(5)

# Merge county data onto shapefile
county_df = pd.merge(county_df,countyPops[['fips','CTYNAME','CTYNAME','POPESTIMATE2019']],on='fips',how='left')

# make population density numbers
county_df['popDensity'] = county_df.POPESTIMATE2019 / county_df.ALAND





######
# Begin analysis
# The example is for Oshkosh, WI on 17 Aug 20
#####

ratioList = [] # collect results of stats
# overlappingList = []
# nationalCloseList = []

# List counties to be excluded (e.g. geographic issues like catching counties across a lake)
excludeList = ['26109','26165','26043','26101','26121', '26123','26085', '26139','26019', '26105', '26127'] # removing the finger of "local" across the great lakesin Michigan

## Commenting out removal of counties "local" to other rallies for simplicity
# for i in allLocalCounties.localRallyCounties.tolist():
#     excludeList.append(i)

# # for row in locs.itertuples():
# for row in locs.itertuples():    
rallyLocation = 'wittman regional airport, oshkosh, wi' #row.address
rallyDay = dt.datetime(2020,8,17,0,0)   #row.rallyDate
location = 'oshkosh_wi' #row.filename
stateOfInterest='Wisconsin'
surroundingCountiesRadius = 150 # number of miles a county's centroid can be from rally location to be included in "surrounding counties"



eventDate = rallyDay.strftime('%-m/%-d/%y') # day of rally

preDate = (rallyDay-dt.timedelta(days=14)).strftime('%-m/%-d/%y') # 14 days before the rally

postDate = (rallyDay+dt.timedelta(days=14)).strftime('%-m/%-d/%y') # 14 days after the rally


# LOCALIZED INFECTION TREND FRAME

# Geocode location
import googlemaps
from datetime import datetime
gmaps = googlemaps.Client(key=googleMapsKey)  ## input your googlemaps API key

locationsList = [rallyLocation]

prevResultsDict = []

resultsDict = []

googleCount = 0

count = 0
for i in range(0,len(locationsList)):
    loc = locationsList[i]

    geocode_result = gmaps.geocode(loc)
    googleCount = googleCount + 1
    locDict = {'inputAddress':loc,'googleResult':geocode_result}
    resultsDict.append(locDict)

    count = count + 1
try:
    c = pd.DataFrame(prevResultsDict)
except:
    c = pd.DataFrame()

r = pd.DataFrame(resultsDict)

allLocs = []

for i in range(0,len(r)):
    inputAddress = r.inputAddress[i]
    if len(r.googleResult[i])>0:
        googleMatchedAddress = r.googleResult[i][0].get('formatted_address')
        lat = r.googleResult[i][0].get('geometry').get('location').get('lat')
        long = r.googleResult[i][0].get('geometry').get('location').get('lng')
    else:
        googleMatchedAddress = 'noMatch'
        lat = 'noMatch'
        long = 'noMatch'
    locDict = {'inputAddress':inputAddress,'googleMatchedAddress':googleMatchedAddress,'lat':lat,'long':long}
    allLocs.append(locDict)
dg = pd.DataFrame(allLocs)
frames = [c,dg]
dg = pd.concat(frames)

# get the state of the input address from the google result
for i in geocode_result[0]['address_components']:
    if i['types'][0]=='administrative_area_level_1':
        stateOfAddress = i['long_name']
lat = dg[0:1].lat.item()  # rally location lattitude
long = dg[0:1].long.item() # rally location longitude


# Use OpenRouteService API to get the 1-hr drive time isochrone
import requests
try:
    body = {"locations":[[long,lat]],"range":[3600],"smoothing":5} # 3600 seconds = 1hr

    headers = {
        'Accept': 'application/json, application/geo+json, application/gpx+xml, img/png; charset=utf-8',
        'Authorization': openrouteserviceKey,    ## Input your Openroutesuervice API key here
        'Content-Type': 'application/json; charset=utf-8'
    }
    call = requests.post('https://api.openrouteservice.org/v2/isochrones/driving-car', json=body, headers=headers)


    iso = json.loads(call.text)
    loc_dict = {}
    name = 'venue'
    fillColor= 'blue'
    loc_dict[name]= {'location':[long,lat],fillColor:fillColor}
    polys = []
    for key in loc_dict.keys():
        geom = shapely.geometry.shape(iso['features'][0]['geometry'])
        polys.append(geom)

    polygons = polys
    boundary = gpd.GeoSeries(cascaded_union(polygons))
except:
    print('openrouteservice failed')
#     # if openrouteservice fails, go to mapbox
#     firstPart = 'https://api.mapbox.com/isochrone/v1/mapbox/driving/'
#     lastPart =  # input API key
#     url = firstPart + str(long)+','+str(lat) + lastPart
#     call = requests.get(url)


#     iso = json.loads(call.text)
#     loc_dict = {}
#     name = 'venue'
#     fillColor= 'blue'
#     loc_dict[name]= {'location':[long,lat],fillColor:fillColor}
#     polys = []
#     for key in loc_dict.keys():
#         geom = shapely.geometry.shape(iso['features'][0]['geometry'])
#         polys.append(geom)

#     polygons = polys
#     boundary = gpd.GeoSeries(cascaded_union(polygons))

# Get counties that overlap with 1-hr isochrone
overlapping = []
for i in range(0,len(county_df)):
    tgt = county_df[i:i+1]
    overlap = shape(tgt.geometry.item()).intersects(shape(boundary.item()))
    if overlap  == True:
        overlapping.append(tgt.GEOID.item())

#
oneHourCounties = county_df[county_df.GEOID.isin(overlapping)]
visualizing_counties = oneHourCounties

c = iso['features'][0]['geometry']['coordinates'][0]
coords = []
for i in c:
    a = [i[1],i[0]]
    coords.append(a)



## Visualize


# function to calculate centroid of a polygon
def getXY(pt):
    return (pt.x, pt.y)
centroidseries = county_df['geometry'].centroid
x,y = [list(t) for t in zip(*map(getXY, centroidseries))]
county_df['x_cent'] = x
county_df['y_cent'] = y

# Make centroids for each county
def distFromLocation(lat,long,y_cent,x_cent):
    return geopy.distance.distance((lat,long),(y_cent,x_cent)).miles
county_df['distance'] = county_df.apply(lambda x: distFromLocation(lat,long,x.y_cent, x.x_cent), axis=1)



# closeList = all counties that overlap with the 1-hour-drive-time polygon

closeList = county_df[(county_df['distance']<surroundingCountiesRadius)&(~county_df.GEOID.isin(overlapping))&(~county_df.GEOID.isin(excludeList))].GEOID.tolist()
closeCounties = county_df[(county_df['distance']<surroundingCountiesRadius)&(~county_df.GEOID.isin(overlapping))&(~county_df.GEOID.isin(excludeList))]


# # Separately generated list of all counties in the US that are "local" to any rally
# doubledCounties = county_df[(county_df['distance']<surroundingCountiesRadius)&(~county_df.GEOID.isin(overlapping))&(county_df.GEOID.isin(allLocalCounties.localRallyCounties.tolist()))]


##
# Calculate Risk Ratios: proportion of people in "local to event" population infected in the time window divided by 
#                        proportion of people in "sorrounding counties" infected in the time window 
# Methodology: Principles of Epidemiology in Public Health Practice, Third Edition An Introduction to Applied Epidemiology and Biostatistics (https://www.cdc.gov/csels/dsepd/ss1978/lesson3/section5.html)
##

# dateFrame = [ FIPS, county_population, countyCumulativeTotal-14-days-before-rally, countyCumulativeTotal-day-of-rally,countyCumulativeTotal-14-days-after-rally]
dateFrame = covidData[['FIPS',preDate,eventDate,postDate]]
dateFrame = pd.merge(dateFrame,county_df[['fips','POPESTIMATE2019']],left_on='FIPS',right_on='fips',how='left')

closeDateFrame = dateFrame[dateFrame.FIPS.isin(closeList)]
localDateFrame = dateFrame[dateFrame.FIPS.isin(overlapping)]

# send each to .csv for review
closeDateFrame.to_csv(os.getcwd() + '/surroundingCounties_rrCalcs_' + location + '.csv',index=False)
localDateFrame.to_csv(os.getcwd() + '/localCounties_rrCalcs_' + location + '.csv',index=False)

# Evaluate relative risk in the 14 days BEFORE to the rally
notExposedCases = closeDateFrame[eventDate].sum() - closeDateFrame[preDate].sum() # cumulative cases outside in surrounding counties day of rally - cumulative cases there 14 days prior to rally
notExposedNotCases = closeDateFrame.POPESTIMATE2019.sum() - closeDateFrame[eventDate].sum() # total population in surrounding counties minus their cumulative infection total on day of rally
exposedCases=localDateFrame[eventDate].sum()-localDateFrame[preDate].sum() # cumulative cases in local counties on day of rally minus cumulative cases in local counties 14 days prior
exposedNotCases = localDateFrame.POPESTIMATE2019.sum()-localDateFrame[eventDate].sum() # total population in local counties (1-hr polygon overlap) minus cumulative total infections for those counties on day of event
try:
    rrPre = (exposedCases / (exposedCases + exposedNotCases)) / (notExposedCases /(notExposedCases + notExposedNotCases))

except:
    rrPre=0

# Evaluate relative risk ratio in the 14 days AFTER the rally
notExposedCases = closeDateFrame[postDate].sum() - closeDateFrame[eventDate].sum()
notExposedNotCases = closeDateFrame.POPESTIMATE2019.sum() - closeDateFrame[postDate].sum()
exposedCases=localDateFrame[postDate].sum()-localDateFrame[eventDate].sum()
exposedNotCases = localDateFrame.POPESTIMATE2019.sum()-localDateFrame[postDate].sum()
try:
    rr = (exposedCases / (exposedCases + exposedNotCases)) / (notExposedCases /(notExposedCases + notExposedNotCases))

except:
    rr=0

    
    
    
# evaluate "local counties" population density divided by population density of "surrounding counties"
localDensity = county_df[county_df.fips.isin(overlapping)].POPESTIMATE2019.sum() / county_df[county_df.fips.isin(overlapping)].ALAND.sum()
surroundingDensity = county_df[county_df.fips.isin(closeList)].POPESTIMATE2019.sum() / county_df[county_df.fips.isin(closeList)].ALAND.sum()
relativeDensity = localDensity/surroundingDensity # Population density of "local" compared to "surrounding" counties

entry= {
    'location':rallyLocation,
    'date':eventDate,
    'close_county_count':len(closeList),
    'riskRatioPrior':rrPre,
    'riskRatioPost':rr,
    'exposedCases':exposedCases,
    'exposedNotCases':exposedNotCases,
    'notExposedCase':notExposedCases,
    'notExposedNotCase':notExposedNotCases,
    'relativeDensity':relativeDensity

}
ratioList.append(entry)

# Send results to .csv for review
pd.DataFrame(ratioList).to_csv(os.getcwd() + '/' + location + '_summaryData.csv',index=False)

##
# Create .html map visualization file with Folium
##

# Modified Josh Rec colors (teal and orange)
isoColor = '#3D9092'
countyColor = '#FF9901'

visualizing_counties = oneHourCounties

# Make a folium map HTML file
m = folium.Map(
    location=[lat,long],
    tiles='CartoDB Positron',
    zoom_start=8
)
if stateOfAddress == 'Louisiana':
    countyType = 'parish'
else:
    countyType = 'county'


# Add 1-hr counties to map
for i in range(0,len(visualizing_counties)):


    try:
        t = list(zip(*visualizing_counties[i:i+1].geometry.item().exterior.coords.xy))
        county_coords = []
        for i in t:
            a = [i[1],i[0]]
            county_coords.append(a)

        folium.vector_layers.PolyLine(county_coords, popup=None, color=countyColor,fill='grey').add_to(m)
    except:
        for j in visualizing_counties[i:i+1].geometry.item():
            t = list(zip(*j.exterior.coords.xy))
            county_coords = []
            for i in t:
                a = [i[1],i[0]]
                county_coords.append(a)

            folium.vector_layers.PolyLine(county_coords, popup=None, color=countyColor,fill='grey').add_to(m)
# Draw the 1-hour driving isochrone
folium.vector_layers.PolyLine(coords, popup=None, tooltip="1-hr Drive Radius",color=isoColor).add_to(m)


# add surrounding counties to map
for i in range(0,len(closeCounties)):
    try:
        t = list(zip(*closeCounties[i:i+1].geometry.item().exterior.coords.xy))
        county_coords = []
        for i in t:
            a = [i[1],i[0]]
            county_coords.append(a)

        folium.vector_layers.PolyLine(county_coords, popup=None, color='#5c9daa',fill='grey').add_to(m)
    except:
        for j in closeCounties[i:i+1].geometry.item():
            t = list(zip(*j.exterior.coords.xy))
            county_coords = []
            for i in t:
                a = [i[1],i[0]]
                county_coords.append(a)

            folium.vector_layers.PolyLine(county_coords, popup=None, color='#5c9daa',fill='grey').add_to(m)

folium.vector_layers.PolyLine(coords, popup=None,color=isoColor).add_to(m)




######################### NOTE: REMOVED THIS FOR SIMPLICITY...CAN BE ADDED BACK IN TO EXLCUDE COUNTIES "LOCAL" TO OTHER RALLIES

# # add counties that would be "surroudning counties" but are exclued becasue they are local to anothe rally
# for i in range(0,len(doubledCounties)):
#     try:
#         t = list(zip(*doubledCounties[i:i+1].geometry.item().exterior.coords.xy))
#         county_coords = []
#         for i in t:
#             a = [i[1],i[0]]
#             county_coords.append(a)

#         folium.vector_layers.PolyLine(county_coords, popup=None, color='darkgrey',fill='grey',opacity=0.4).add_to(m)
#     except:
#         for j in doubledCounties[i:i+1].geometry.item():
#             t = list(zip(*j.exterior.coords.xy))
#             county_coords = []
#             for i in t:
#                 a = [i[1],i[0]]
#                 county_coords.append(a)

#             folium.vector_layers.PolyLine(county_coords, popup=None, color='darkgrey',fill='grey',opacity=.4).add_to(m)




# add rally location to map
folium.CircleMarker(
    location=[lat, long],
    radius=5,
    tooltip='Rally location',
    color=isoColor,
    fill=True,
    fill_color='#6bbbed'
).add_to(m)

### If re-adding control for counties that are "local" to other rallies, here is the legend code:     <li><span style='background:""" + 'darkgrey' + """;opacity:0.6; border: 1px solid """ + '#5c9daa' + """;'></span>Counties near another rally</li>


# Add legend to map
from branca.element import Template, MacroElement

template = """
{% macro html(this, kwargs) %}

<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">


</head>
<body>


<div id='maplegend' class='maplegend'
    style='position: absolute; z-index:9999; border:2px solid grey; background-color:rgba(255, 255, 255, 0.8);
     border-radius:6px; padding: 10px; font-size:14px; right: 20px; top: 20px;'>

<div class='legend-scale'>
  <ul class='legend-labels'>
    <li><span style='background:""" + countyColor + """;opacity:0.6; border: 1px solid """ + countyColor + """;'></span>1-hour Counties</li>
    <li><span style='background:""" + '#5c9daa' + """;opacity:0.6; border: 1px solid """ + '#5c9daa' + """;'></span>Close Counties</li>

  </ul>
</div>
</div>

</body>
</html>

<style type='text/css'>
  .maplegend .legend-title {
    text-align: left;
    margin-top: 5px;
    font-weight: bold;
    font-size: 90%;
    }
  .maplegend .legend-scale ul {
    margin: 0;
    margin-top: 5px;
    padding: 0;
    float: left;
    list-style: none;
    }
  .maplegend .legend-scale ul li {
    font-size: 80%;
    list-style: none;
    margin-left: 0;
    line-height: 18px;
    margin-top: 2px;
    }
  .maplegend ul.legend-labels li span {
    display: block;
    float: left;
    height: 16px;
    width: 30px;
    margin-right: 5px;
    margin-left: 0;
    border: 1px solid #999;
    }
  .maplegend .legend-source {
    font-size: 80%;
    color: #777;
    clear: both;
    }
  .maplegend a {
    color: #777;
    }
</style>
{% endmacro %}"""

macro = MacroElement()
macro._template = Template(template)

m.get_root().add_child(macro)


# save map to "results" folder
m.save(os.getcwd()+'/' + location + '_map.html')



###
# Visualize trends in infections
###

## Inputs:
# Rally datetime: rallyDay
# Rally address: rallyLocation
# covidData dataframe
# covidDeaths dataframe

# # NATIONAL INFECTION TREND FRAME

# rallyLocation = "Sturgis, SD" #locs[4:5].address.item()
# rallyDay = dt.datetime(2020,8,7,0,0)#locs[4:5].day.item()
# location = rallyLocation

# def makeDifInDif(rallyLocation, rallyDay,filename):
    


# Make frame # make seven day average dataframe ( 'date': 'average new infections in the 7 days prior') for national
#  rolling 7 day averages nationally
sevenDayAvg = []
for i in range(8,len(dates)):
    d = dates[i]
    wBefore = dates[i-7]
    weeklyNew = covidData[d].sum() - covidData[wBefore].sum()
    weeklyDeaths = covidDeaths[d].sum() - covidDeaths[wBefore].sum()

    #     perHundK = round(100000*dailyNew / countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
    #     activePercenage = round(100*(confirmed[confirmed.FIPS.isin(overlapping)][d].sum() )/countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
    entry = {'date':d,
             'sevenDayAvg':weeklyNew/7,
             'sevenDayAvgDeaths':weeklyDeaths/7}
#             'perHundK':perHundK}
    sevenDayAvg.append(entry)
sevenDayAvg = pd.DataFrame(sevenDayAvg)
sevenDayAvg['day'] = pd.to_datetime(sevenDayAvg.date)

# Cut natonal sevenDay average down to 60 days before and after the date of interest
nationalWindow = sevenDayAvg[(sevenDayAvg.day >= (rallyDay-dt.timedelta(days=60)))&\
           (sevenDayAvg.day <= (rallyDay+dt.timedelta(days=60)))]
# Make a new column of infections normalized to the value on the day of interest
nationalWindow['infectionTrend'] =  100*nationalWindow['sevenDayAvg']/(nationalWindow[nationalWindow.day==rallyDay].sevenDayAvg.item())



# make seven day average dataframe ( 'date': 'average new infections in the 7 days prior') for local counties

localSevenDay = []
for i in range(8,len(dates)):
    d = dates[i]
    wBefore = dates[i-7]
    weeklyNew = covidData[covidData.FIPS.isin(overlapping)][d].sum() - covidData[covidData.FIPS.isin(overlapping)][wBefore].sum()
    weeklyDeaths = covidDeaths[covidDeaths.FIPS.isin(overlapping)][d].sum() - covidDeaths[covidDeaths.FIPS.isin(overlapping)][wBefore].sum()


    #     perHundK = round(100000*dailyNew / countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
    #     activePercenage = round(100*(confirmed[confirmed.FIPS.isin(overlapping)][d].sum() )/countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
    entry = {'date':d,
             'sevenDayAvg':weeklyNew/7,
             'sevenDayAvgDeaths':weeklyDeaths/7}
#             'perHundK':perHundK}
    localSevenDay.append(entry)
localSevenDay = pd.DataFrame(localSevenDay)
localSevenDay['day'] = pd.to_datetime(localSevenDay.date)
# Cut local sevenDay average down to 60 days before and after the date of interest
localWindow = localSevenDay[(localSevenDay.day >= (rallyDay-dt.timedelta(days=60)))&\
           (localSevenDay.day <= (rallyDay+dt.timedelta(days=60)))]
# Make a new column of infections normalized to the value on the day of interest
localWindow['infectionTrend'] =  100*localWindow['sevenDayAvg']/(localWindow[localWindow.day==rallyDay].sevenDayAvg.item())

## Visualize



# make seven day average dataframe ( 'date': 'average new infections in the 7 days prior') for surrounding counties
closeSevenDayAvg = []
for i in range(8,len(dates)):
    d = dates[i]
    wBefore = dates[i-7]
    weeklyNew = covidData[covidData.FIPS.isin(closeList)][d].sum() - covidData[covidData.FIPS.isin(closeList)][wBefore].sum()
    weeklyDeaths = covidDeaths[covidData.FIPS.isin(closeList)][d].sum() - covidDeaths[covidData.FIPS.isin(closeList)][wBefore].sum()

    #     perHundK = round(100000*dailyNew / countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
    #     activePercenage = round(100*(confirmed[confirmed.FIPS.isin(overlapping)][d].sum() )/countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
    entry = {'date':d,
             'sevenDayAvg':weeklyNew/7,
             'sevenDayAvgDeaths':weeklyDeaths/7}
#             'perHundK':perHundK}
    closeSevenDayAvg.append(entry)
closeSevenDayAvg = pd.DataFrame(closeSevenDayAvg)
closeSevenDayAvg['day'] = pd.to_datetime(closeSevenDayAvg.date)
closeSevenDayAvg['infectionTrend'] =  100*closeSevenDayAvg['sevenDayAvg']/(closeSevenDayAvg[closeSevenDayAvg.day==rallyDay].sevenDayAvg.item())

closeWindow = closeSevenDayAvg[(closeSevenDayAvg.day >= (rallyDay-dt.timedelta(days=60)))&\
           (closeSevenDayAvg.day <= (rallyDay+dt.timedelta(days=60)))]




# ## Make state rolling-7-day-average trend (commented out for visual clarity)
# stateSevenDayAvg = []
# for i in range(8,len(dates)):
#     d = dates[i]
#     wBefore = dates[i-7]
#     weeklyNew = covidData[covidData.Province_State==stateOfInterest][d].sum() - covidData[covidData.Province_State==stateOfInterest][wBefore].sum()
#     weeklyDeaths = covidDeaths[covidDeaths.Province_State==stateOfInterest][d].sum() - covidDeaths[covidDeaths.Province_State==stateOfInterest][wBefore].sum()

#     #     perHundK = round(100000*dailyNew / countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
#     #     activePercenage = round(100*(confirmed[confirmed.FIPS.isin(overlapping)][d].sum() )/countyPops[countyPops.FIPS.isin(overlapping)].POPESTIMATE2019.sum(),2)
#     entry = {'date':d,
#              'sevenDayAvg':weeklyNew/7,
#              'sevenDayAvgDeaths':weeklyDeaths/7}
# #             'perHundK':perHundK}
#     stateSevenDayAvg.append(entry)
# stateSevenDayAvg = pd.DataFrame(stateSevenDayAvg)
# stateSevenDayAvg['day'] = pd.to_datetime(stateSevenDayAvg.date)

# stateWindow = stateSevenDayAvg[(stateSevenDayAvg.day >= (rallyDay-dt.timedelta(days=60)))&\
#            (stateSevenDayAvg.day <= (rallyDay+dt.timedelta(days=60)))]
# stateWindow['infectionTrend'] =  100*stateWindow['sevenDayAvg']/(stateWindow[stateWindow.day==rallyDay].sevenDayAvg.item())


plt.figure(figsize=(20,10),dpi=300)
ax = plt.subplot(111)
# plt.figsize=[.2,.4]



plt.fill_between(nationalWindow.day.tolist(), nationalWindow.infectionTrend.tolist(), 100, facecolor='lightgrey', alpha=0.5)

# daily.set_index('dt').perHundK.plot.line(linewidth=4,color=countyColor)
localWindow.set_index('day').infectionTrend.plot.line(linewidth=7,color=countyColor)
nationalWindow.set_index('day').infectionTrend.plot.line(linewidth=4,color='grey')
closeWindow.set_index('day').infectionTrend.plot.line(linewidth=4,color='#5c9daa')
# stateWindow.set_index('day').infectionTrend.plot.line(linewidth=4,color='red')

# plt.title("Percentage of Population Infected \nwithin an Hour's Drive of {}\n\n".format(location),fontsize=25)
if len(location)>15:
    plt.title(location.split(',')[0] + '\n' + location.split(',')[1] + ', ' + location.split(',')[2] + '\n \n',fontsize=50)
else:
    plt.title(location + '\n \n',fontsize=50)
#     figtext(.5,.9,'Foo Bar', fontsize=18, ha='center')

plt.xlabel('')

ax.set_ylim([-50,700])

# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
plt.yticks([])

# ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')

# percInfectedList = weekly.perHundK.tolist()




# ax.tick_params(axis="y", labelsize=15,color='grey')
ax.tick_params(axis="x", labelsize=20,color='grey')

plt.minorticks_off()

# plt.text(max(daily.dt+dt.timedelta(days=2)),daily.perHundK.max(),maxText,horizontalalignment='left',fontsize=20,color='grey')

# Work end label for local infection trend

if localWindow[-1:].infectionTrend.item()-100>0:
#         lastDay = '{:,}'.format(round(localWindow[-1:].infectionTrend.item()-100,0))
    lastDay = str(int(round(localWindow[-1:].infectionTrend.item()-100,0)))

    labelText = 'Local infections up ' +lastDay +'%'
    labelTwoText =  str(localWindow[-1:].infectionTrend.item()) + ' % up'

    plt.text(max(localWindow.day+dt.timedelta(days=1)),localWindow.infectionTrend[-1:].item(),labelText,horizontalalignment='left',fontsize=30)

else:
    downNumber = int(100 - round(localWindow[-1:].infectionTrend.item()))
    lastDay = str(downNumber)

    labelText = 'Local infections down ' +lastDay +'%'
#         labelTwoText =  str(localWindow[-1:].infectionTrend.item()) + ' % up'

    plt.text(max(localWindow.day+dt.timedelta(days=1)),localWindow.infectionTrend[-1:].item(),labelText,horizontalalignment='left',fontsize=30)

# Work end label for national infection trend

if nationalWindow[-1:].infectionTrend.item()-100>0:
    lastDay = str(int(round(nationalWindow[-1:].infectionTrend.item()-100,0)))

    labelText = 'National infections up ' +lastDay +'%'
#         labelTwoText =  str(localWindow[-1:].infectionTrend.item()) + ' % up'

    plt.text(max(nationalWindow.day+dt.timedelta(days=1)),nationalWindow.infectionTrend[-1:].item(),labelText,horizontalalignment='left',fontsize=30)

else:
    downNumber = int(100 - round(nationalWindow[-1:].infectionTrend.item()))
    lastDay = str(downNumber)

    labelText = 'National infections down ' +lastDay +'%'
#         labelTwoText =  str(localWindow[-1:].infectionTrend.item()) + ' % up'

    plt.text(max(nationalWindow.day+dt.timedelta(days=1)),nationalWindow.infectionTrend[-1:].item(),labelText,horizontalalignment='left',fontsize=30)


if closeWindow[-1:].infectionTrend.item()-100>0:
    lastDay = str(int(round(closeWindow[-1:].infectionTrend.item()-100,0)))

    labelText = 'Surrounding counties infections up ' +lastDay +'%'
#         labelTwoText =  str(localWindow[-1:].infectionTrend.item()) + ' % up'

    plt.text(max(closeWindow.day+dt.timedelta(days=1)),closeWindow.infectionTrend[-1:].item(),labelText,horizontalalignment='left',fontsize=30)
else:
    downNumber = int(100 - round(nationalWindow[-1:].infectionTrend.item()))
    lastDay = str(downNumber)

    labelText = 'Surrounding counties infections down ' +lastDay +'%'

    

## Add state trend (commented out for clarity of visual presentation)
# if stateWindow[-1:].infectionTrend.item()-100>0:
#     lastDay = str(int(round(stateWindow[-1:].infectionTrend.item()-100,0)))

#     labelText = 'State infections up ' +lastDay +'%'
# #         labelTwoText =  str(localWindow[-1:].infectionTrend.item()) + ' % up'

#     plt.text(max(stateWindow.day+dt.timedelta(days=1)),stateWindow.infectionTrend[-1:].item(),labelText,horizontalalignment='left',fontsize=30)
# else:
#     downNumber = int(100 - round(stateWindow[-1:].infectionTrend.item()))
#     lastDay = str(downNumber)

#     labelText = 'State infections down ' +lastDay +'%'    
    
    
# add label at top of "event date" line
plt.text(rallyDay-dt.timedelta(days=15),680,rallyDay.strftime('%d %b %Y') + '\n',fontsize=35,color='grey')                       



plt.axvline(x=rallyDay, color='k', linestyle='--')

plt.axhline(y=100, color='lightgrey', linestyle='--')


plt.tight_layout()

fn = os.getcwd() +'/' + location + '_trend.png'
plt.savefig(fn, dpi=300, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        metadata=None)

