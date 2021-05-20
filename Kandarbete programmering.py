# -*- coding: utf-8 -*-


import cartopy.crs as ccrs
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import basemap
from netCDF4 import Dataset
import netCDF4
from matplotlib.patches import Rectangle
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
from scipy import stats

#Datum: 2021-05-16
#Kod för kandidatarbete: Stratocumulus på låga latituder: variabilitet och relevans för klimatet
#Författare: Onni Mikkola

"""
Kodbeskrivning: 
Översta delen av koden består av funktioner som programmet kallar på. 
Programmet som kallar på funktionerna ligger längst ner i koden. 


Output från programmet: 
- Global fördelning av molnfraktionen av alla moln och alla låga moln.
- Olika spridningsdiagram mellan molnfraktion och stabilitetsmåttet LTS.
- Tidsutveckling av LTS


Vad som behövs för att kod skall fungera:
Programtolk: Python 2.7
Moduler (kan importeras via anaconda t.ex.): matplotlib, numpy, cartopy, mpl_toolkits, netCDF4
Nedladdade filer skall ligga i samma map som denna kod.


ERA5 datan har laddats ner från https://cds.climate.copernicus.eu/cdsapp#!/search?type=dataset [2021-05-16]
där följande data för laddats ner:
- ERA5 monthly averaged data on single levels from 1979 to present
- ERA5 monthly averaged data on pressure levels from 1979 to present
#Datan har laddats ner för de beskrivna latituderna/longituderna i rapporten
    - Den nedladdade datan är temperaturen på 1000 och 700 hPa samt molnfraktionen av låga moln.
#Samt en global fördelning för 2017-2020
    - Den nedladdade datan är molnfraktionen av låga moln.

För att koden skall fungera bör filerna döpas till:
Pressurelvl_79-20_Peruvian (Samma för: Namibian, Australian, Californian, Canarian)
Singlelvl_79-20_Peruvian (Samma för: Namibian, Australian, Californian, Canarian)
Singlelvl_17-20_Global
"""


#Creates map projection for plot
def cartopy_map(projection):
    # This function calls the desired projection for the cartopy.

    #Flat projection
    if projection == ccrs.PlateCarree:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()
        ax.gridlines()
        ax.set_xticks([-180,-120 ,-60 ,0, 60, 120, 180], crs=ccrs.PlateCarree())
        ax.set_yticks([-80 , -40 , 0 , 40 ,80], crs=ccrs.PlateCarree())
        ax.set_xticklabels([-180,-120 ,-60 ,0, 60, 120, 180], fontsize=30)
        ax.set_yticklabels([-80 , -40 , 0 , 40 ,80], fontsize=30)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

    #Globe projection
    if projection == ccrs.Orthographic:
        # Add Orthographic projection here.
        ax = plt.axes(projection=ccrs.Orthographic(central_longitude=5, central_latitude=-15, globe=None))
        ax.coastlines()
        ax.set_global()
        ax.gridlines()

    return ax


#Creates global distribution plots
def makePlot(projection, field, lat, lon, fig_title, title, levels, cmap):
    # This function makes a lat/lon plot of given field.
    ###########################################################################
    # Projection = The projection of basemap you want. Used in base_map function.
    # field      = The field you want to plot.
    # lat        = latitude of the field.
    # lon        = longitude of the field.
    # fig_title  = title of the figure.
    # levels     = putting the values of your plotting limits.
    # cmap       = The colormap you would like to use for the specific plot.
    ###########################################################################
    fig = plt.figure(fig_title)
    w = 16
    h = 8
    fig.set_size_inches(w, h)
    fig.patch.set_facecolor('white')
    m = cartopy_map(projection)  # This calls the previous function to
    # generate the projection of choosen cartopy map.
    field, lon = basemap.addcyclic(field, lon)
    field, lon = basemap.shiftgrid(180, field, lon, start=False)
    lon, lat = np.meshgrid(lon, lat)
    cf = m.contourf(lon, lat, field, 20, cmap=cmap, levels=levels, extend='both', transform=ccrs.PlateCarree())
    cb = plt.colorbar(cf, ticks=levels[::2])
    cb.ax.tick_params(labelsize=30)
    currentAxis = plt.gca()
    regionleftcorners = [0,-20,-90,-20,-130,20,95,-35,-35,15]
    for i in range(0,len(regionleftcorners),2):
        currentAxis.add_patch(Rectangle((regionleftcorners[i],regionleftcorners[i+1],), 10, 10,fc = "None", color ='black',linewidth = 2,linestyle="-"))
    plt.title(title, fontsize = 32)
    plt.tight_layout(pad=0.3, h_pad=0.3, w_pad=0.3)
    plt.show()


#Read nc file
def readfilepslevel(filename):
    filepath = "C:/Users/Onni/Desktop/Koulu/Kandidatprogram i Meteorologi/Arskurs 3/Bachelor thesis/Programmering/" + filename + ".nc"
    datapslevel = Dataset(filepath, 'r')
    return datapslevel


#Calculates LTS for a region
def CALCLTS(datapslevel,regionname, annual):

    #Extract temperature at 1000 & 700 hPa
    T1000 = datapslevel.variables['t'][:, 1, :, :]  # 3 axes: time,lat,lon
    T700 = datapslevel.variables['t'][:, 0, :, :]

    if annual == True:

        # Calculates PotT at 1000 & 700 hPa [C] (Annual mean for whole timeseries)
        PotT1000 = T1000.mean(0).mean(0).mean(0) - 273.15
        PotT700 = (T700.mean(0).mean(0).mean(0) * (float(1000) / 700) ** 0.286) - 273.15

        # Calculates LTS in each grid
        LTS = PotT700 - PotT1000

        print("The annual mean LTS of " + regionname + " is: " + str(LTS))

        return LTS, T1000, T700, PotT1000, PotT700

    elif annual == False:
        allaar1000 = []
        allaar700 = []
        LTSmonth = []

        # go through every unique month
        for i in range(0, 12):

            #Creates empty list for every month
            uniktarslista1000 = []
            uniktarslista700 = []

            # k = what month of the year
            k = i

            #Go through every year, len(T1000)/12 = amount of years.
            for j in range(0, len(T1000)/12):

                #Add January e.g. to list for all years
                uniktarslista1000.append(T1000[k])
                uniktarslista700.append(T700[k])

                #Go 12 steps forward to get to next years january
                k += 12

            allaar1000.append(uniktarslista1000)
            allaar700.append(uniktarslista700)

        allaar1000 = np.array(allaar1000)
        allaar700 = np.array(allaar700)

        for l in range(0,12):
            allaar1000m = allaar1000[l].mean(0).mean(0).mean(0)
            allaar700m = allaar700[l].mean(0).mean(0).mean(0)

            #Calculates the mean potential temp for whole area per month
            PotT1000 = allaar1000m-273.15
            PotT700 = (allaar700m * (float(1000)/700)**0.286) - 273.15

            #Calculates mean LTS per area per month
            LTS = PotT700 - PotT1000

            LTSmonth.append(LTS)

            monthlist = ['january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september','october', 'november', 'december']
            print("The mean LTS in " + monthlist[l] + " in " + regionname + " is: " + str(LTS))

        k = -1
        lista1000 = []
        lista700 = []

        #Calculates LTS for djf, mam...
        for i in range(0,4):

            monthlymeans1000 = (1/float(3))*(allaar1000[k].mean(0).mean(0).mean(0) + allaar1000[k+1].mean(0).mean(0).mean(0) + allaar1000[k+2].mean(0).mean(0).mean(0))
            monthlymeans700= (1/float(3))*(allaar700[k].mean(0).mean(0).mean(0) + allaar700[k+1].mean(0).mean(0).mean(0) + allaar700[k+2].mean(0).mean(0).mean(0))
            lista1000.append(monthlymeans1000)
            lista700.append(monthlymeans700)
            k += 3

        lista1000 = np.array(lista1000)
        lista700 = np.array(lista700)

        LTSnotsat = []

        for i in range(0,4):

            # Calculates the mean potential temp for whole area per month
            PotT1000 = lista1000[i] - 273.15
            PotT700 = (lista700[i] * (float(1000) / 700) ** 0.286)-273.15

            # Calculates mean LTS per area per month
            LTS = PotT700 - PotT1000

            LTSnotsat.append(LTS)

            months = ['djf', 'mam', 'jja', 'son']
            print("The mean LTS in " + str(months[i]) + " in " + regionname + " is: " + str(LTS))


    global lista1000, lista700, PotT1000, PotT700,allaar1000,allaar700, LTSmonth, T1000, T700, LTSnotsat
    return LTS, allaar1000,allaar700 ,PotT1000,PotT700, LTSmonth, T1000, T700


#Calculates CF (cloud fraction) for a region
def CALCCF(datapslevel,regionname, annual):

    # Extract variables from nc file
    LCC = datapslevel.variables['lcc'][:, :, :]  # 3 axes: time,lat,lon
    # TCC = datapslevel.variables['tcc'][:, :, :]  # 3 axes: time,lat,lon

    CFnotsat = []

    if annual == True:

        #Calculate LCC mean
        LCCm = LCC.mean(0).mean(0).mean(0)

        print("The annual mean CF of " + regionname + " is: " + str(LCCm))

    elif annual == False:
        allaarCF = []
        allaarCFmonth = []

        # Go through every unique month
        for i in range(0, 12):

            #Create empty list for every month
            uniktarslistaCF = []

            # k = what month of the year
            k = i

            #Go through every year,len(LCC)/12 = amount of years.
            for j in range(0, len(LCC)/12):
                # lagg till varje januari exempelvis till listan for alla ar
                #Add every january e.g. to list for all years
                uniktarslistaCF.append(LCC[k])

                #Go 12 steps forward to get to next january
                k += 12

            allaarCF.append(uniktarslistaCF)

        allaarCF = np.array(allaarCF)


        for l in range(0,12):
            allaarCFm = allaarCF[l].mean(0).mean(0).mean(0)

            allaarCFmonth.append(allaarCFm)

            monthlist = ['january','february','march','april','may','june','july','august','september','october','november','december']
            print("The mean CF in " + monthlist[l] + " in " + regionname + " is: " + str(allaarCFm))

        k = -1
        listaCF = []

        #Calculates cloud fraction for djf, mam...
        for i in range(0,4):
            monthlymeansCF = (1/float(3))*(allaarCF[k].mean(0).mean(0).mean(0) + allaarCF[k+1].mean(0).mean(0).mean(0) + allaarCF[k+2].mean(0).mean(0).mean(0))
            listaCF.append(monthlymeansCF)
            k += 3

        CFnotsat.append(listaCF)
        listaCF = np.array(listaCF)

        for i in range(0,4):
            months = ['djf', 'mam', 'jja', 'son']
            print("The mean CF in " + str(months[i]) + " in " + regionname + " is: " + str(listaCF[i]))

    global allaarCFm, allaarCF, listaCF, allaarCFmonth, LCC, CFnotsat
    return LCC


#Function that decides what function to call depending on input argument
def LTSCF(filename,regionname, LTS, annual):
    if annual == True:

        if LTS == True:
            datapslevel = readfilepslevel(filename)
            CALCLTS(datapslevel,regionname, annual)
            return

        elif LTS == False:
            datapslevel = readfilepslevel(filename)
            CALCCF(datapslevel, regionname, annual)
            return

    elif annual == False:
        if LTS == True:
            datapslevel = readfilepslevel(filename)
            CALCLTS(datapslevel,regionname, annual)
            return

        elif LTS == False:
            datapslevel = readfilepslevel(filename)
            CALCCF(datapslevel, regionname, annual)
            return


#Creates global distribution of low- and total cloud amount
def globaldistrib(filenameglobaldistrib):
    #Read nc file
    datapslevel = readfilepslevel(filenameglobaldistrib)

    #Extract variables from nc file
    LCC = datapslevel.variables['lcc'][:, :, :]  # 3 axes: time,lat,lon
    TCC = datapslevel.variables['tcc'][:, :, :]  # 3 axes: time,lat,lon
    lat = datapslevel.variables['latitude'][:]
    lon = datapslevel.variables['longitude'][:]

    #Calculate LCC/TCC mean
    #LCCm = LCC[:, :, :]  # select all januarys
    LCCm = LCC.mean(0)
    TCCm = TCC.mean(0)


    #Define colorbar intervall and call for plot functions
    levels  = np.linspace(0, 1, 21)  # min, max, spacing
    makePlot(ccrs.PlateCarree, LCCm, lat, lon, u"Årligt medel av molnfraktion av låga moln (2017-2020)", u"Årligt medel av molnfraktion av låga moln (2017-2020)", levels,
         'RdYlBu_r')
    makePlot(ccrs.PlateCarree, TCCm, lat, lon, 'Annual Total Cloud Fraction (2017-2020)', 'Annual Total Cloud Fraction (2017-2020)', levels,
         'RdYlBu_r')


#Creates timeseries plot of LTS to time in desired timeintervall between 1979-2020
def timeseries(regionlist,startyear,endyear):
    startyear = startyear * 12
    endyear = endyear * 12
    regionlista = ['Peru','Namibia', 'Kalifornien','Australien',u'Kanarieöarna']

    # fontsizes
    legendsize = 25
    titlesize = 35
    xtitlesize = 26
    ytitlesize = 26
    xtickssize = 25
    ytickssize = 25
    linesize = 2

    #Iterate each region
    for i in range(len(regionlist)):
        #Read nc file
        datatimeseries = readfilepslevel("Pressurelvl_79-20_" + regionlist[i])

        #Extract variable from nc file
        T1000timeseries = datatimeseries.variables['t'][:, 1, :, :]  # 4 axes: time,pressurelevel,lat,lon
        T700timeseries = datatimeseries.variables['t'][:, 0, :, :]
        ttimeseries = datatimeseries.variables['time']
        dtime = netCDF4.num2date(ttimeseries[startyear:endyear], ttimeseries.units) #convert time into datetime

        #Slice arrays into desired timeperiod
        T1000timeseriessliced = T1000timeseries[startyear:endyear].mean(1).mean(1)
        T700timeseriessliced = T700timeseries[startyear:endyear].mean(1).mean(1)
        ttimeseries = ttimeseries[startyear:endyear]

        #Calculate mean LTS for each month between end- and startyear
        PotT1000timeseries = T1000timeseriessliced - 273.15
        PotT700timeseries = (T700timeseriessliced * (float(1000) / 700) ** 0.286) - 273.15
        LTStimeseries = PotT700timeseries - PotT1000timeseries

        ### Seasonal adjustment###

        #empty lists and definitions for calculations below
        allyearaverages = []
        seasonalindecies = []
        monthseasonalaverages = []
        years = (endyear - startyear)/12
        k = 0

        #Calculates seasonalindex for each datapoint and
        for p in range(0, years):
            allyearaverages.append(np.sum(LTStimeseries[k:k + 12])/float(12)) #average LTS per year

            # Seasonal index per datapoint. Divide datapoint by average LTS that year.
            seasonalindecies.append(LTStimeseries[k:k + 12] / float(allyearaverages[p]))
            k += 12

        #Convert seasonalindecies list into a 1-dim array
        seasonalindecies = np.array(seasonalindecies)
        seasonalindecies = np.hstack(seasonalindecies)

        #Calculate monthly average seasonalindex for whole timeperiod
        for p in range(0, 12):
            monthseasonalaverages.append(np.sum(seasonalindecies[p::12]/years))

        #Divide timeseries data by average seasonal index for that month
        deseasonalizationallyears = []
        k = 0
        for p in range(0, years):
            deseasonalizationayear = []
            for j in range(0, 12):
                deseasonalizationayear.append(LTStimeseries[k + j] / float(monthseasonalaverages[j]))
            k += 12
            deseasonalizationallyears.append(deseasonalizationayear)

        #Convert deseasonalization data into 1-dim array
        deseasonalization = np.array(deseasonalizationallyears)
        deseasonalization = np.hstack(deseasonalization)

        x = ttimeseries
        y = deseasonalization

        def fit_func(x, a, b):
            return a * x + b

        # What the values are:popt[0]=a, popt[1]=b, pcov[0,0]=a_error, pcov[1,1]=b_error
        params = curve_fit(fit_func, x, y)  # Curvefit estimates a and b
        [a, b] = params[0]
        popt, pcov = curve_fit(fit_func, x, y)

        legendlist = ['LTS',u'Säsongsjusterad LTS']
        legendlist.append('Trendlinje')

        """
        legendlist.append(
            'Trendlinje (y = ax + b) \n ' +
            u'LTS ökning/år: \n'
            + " a = " + str(round(popt[0]*365*24, 3)) + " $\pm$ " + str(round((pcov[0, 0] ** 0.5)*365*24, 3))
            +' K')
        """

        #Calculate and print p-value for each region
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, LTStimeseries)
        print('p-value for ' + regionlist[i] + ' is: ' + str(p_value))


        #Creates plot
        plt.figure()
        LTStimeseries1, = plt.plot(dtime,LTStimeseries)
        deseasonalization, = plt.plot(dtime, deseasonalization)
        deseasonalizedLTS, = plt.plot(dtime,a*x+b, color = 'red', linewidth = linesize)
        plt.title(u"Tidsserie för LTS i " + regionlista[i] + ' (' + str(1979 + startyear/12) + '-' + str(1979 + endyear/12) + ')', size = titlesize)
        plt.xlabel(u'År', size = xtitlesize)
        plt.ylabel('LTS [K]', size = ytitlesize)
        plt.xticks(rotation = 15, size = xtickssize)
        plt.yticks(size=ytickssize)
        plt.xlim(min(dtime), max(dtime))
        plt.legend([LTStimeseries1,deseasonalization,deseasonalizedLTS],legendlist, fontsize = legendsize)
        plt.show()

        global LTStimeseries1, LTStimeseries, dtime, deseasonalization, x, y


#Creates scatterplot between different kinds of variables
def scatterplot(LTSmonthvalues,CFmonthvalues, satellitedata):
    x = LTSmonthvalues
    y = CFmonthvalues
    regionnames = ['P', 'N', 'Cal', 'A', 'Ca']
    regionlista = ['Peru','Namibia', 'Kalifornien','Australien',u'Kanarieöarna']

    #fontsizes
    legendsize = 30
    titlesize = 35
    monthsize = 25
    regionsize = 25
    xtitlesize = 30
    ytitlesize = 30
    xtickssize = 30
    ytickssize = 30
    linesize = 4
    datapointsize = 15
    datapointsize1 = 10


    def fit_func(x, a, b):
        return a * x + b

    # What the values are:popt[0]=a, popt[1]=b, pcov[0,0]=a_error, pcov[1,1]=b_error
    params = curve_fit(fit_func, x, y)  # Curvefit estimates a and b
    [a, b] = params[0]
    popt, pcov = curve_fit(fit_func, x, y)

    # Have to take line end points to get bestfit dashed in plot. Datapoints not sorted.
    point1 = [min(x) - 1, a * (min(x) - 1) + b]
    point2 = [max(x) + 1, a * (max(x) + 1) + b]
    linexvalues = [point1[0], point2[0]]
    lineyvalues = [point1[1], point2[1]]

    # calculate Pearson's correlation
    corr, _ = pearsonr(x, y)

    if satellitedata == 'seasonsat' or satellitedata == 'seasonnotsat':
        legendlist = ['DJF','MAM','JJA','SON']

    else:
        legendlist = ['Datapunkter']

    """
    legendlist.append(
        'Trendlinje (y = ax + b) \n'
        + " a = " + str(round(popt[0], 3)) + " $\pm$ " + str(round(pcov[0, 0] ** 0.5, 3)) + "\n"
        #+ " b = " + str(round(popt[1], 2)) + " $\pm$ " + str(round(pcov[1, 1] ** 0.5, 2)) + "\n"
        + " r = %.2f" % corr)
    """

    # Creates scatterplot
    plt.figure()

    ###Creates different kind of scatterplot depending on input argument###
    if satellitedata == True:
        i = h
        satelliteobserved = u' (satellitobserverad målnfraktion)'
        bestfit, = plt.plot(linexvalues, lineyvalues, linestyle = 'dashed',linewidth = linesize, color='grey')
        datapoints, = plt.plot(x, y, 'o', color='red',markersize = datapointsize, markerfacecolor='red')
        plt.title(u'Regressionsanalys för ' + regionlista[i] + satelliteobserved + ',' + u'\n medel av månadsmedel (2001-2020)', size=titlesize)
        plt.ylim(min(y)-0.1, max(y)+0.1)
        #plt.legend([datapoints, bestfit], legendlist, fontsize= legendsize)
        #plt.legend([datapoints], legendlist, fontsize=legendsize)
        for l in range(12):
            plt.text(x[l], y[l], str(monthnames[l]), fontsize=monthsize)

    elif satellitedata == 'allmonthssat':
        i = h
        satelliteobserved = u' (satellitobserverad molnfraktion)'
        bestfit, = plt.plot(linexvalues, lineyvalues, linestyle='dashed',linewidth = linesize, color='grey')
        datapoints, = plt.plot(x, y, 'o', color='red', markersize=datapointsize1 ,markerfacecolor='red')
        plt.title(u'Regressionsanalys för ' + regionlista[i] + satelliteobserved + ',' + u'\n månadsmedel (2001-2020)', size=titlesize)
        plt.ylim(min(y) - 0.1, max(y) + 0.1)
        #plt.legend([datapoints, bestfit], legendlist, fontsize=legendsize)

    elif satellitedata == False:
        i = s
        satelliteobserved = ''
        bestfit, = plt.plot(linexvalues, lineyvalues, linestyle='dashed',linewidth = linesize, color='grey')
        datapoints, = plt.plot(x, y, 'o', color='red', markerfacecolor='red',markersize=datapointsize)
        plt.title(u'Regressionsanalys för ' + regionlista[i] + satelliteobserved + ',' + u'\n medel av månadsmedel (1979-2020)', size=titlesize)
        plt.ylim(min(y) - 0.1, max(y) + 0.1)
        #plt.legend([datapoints, bestfit], legendlist, fontsize=legendsize)
        for l in range(12):
            plt.text(x[l], y[l], str(monthnames[l]), fontsize=monthsize)

    elif satellitedata == 'allmonthsnotsat':
        i = s
        satelliteobserved = ''
        bestfit, = plt.plot(linexvalues, lineyvalues, linestyle='dashed',linewidth = linesize, color='grey')
        datapoints, = plt.plot(x, y, 'o', color='red',markersize = datapointsize1, markerfacecolor='red')
        plt.title(u'Regressionsanalys för ' + regionlista[i] + satelliteobserved + ',' + u'\n månadsmedel (1979-2020)', size=titlesize)
        plt.ylim(min(y) - 0.1, max(y) + 0.1)
        #plt.legend([datapoints, bestfit], legendlist, fontsize=legendsize)

    elif satellitedata == 'seasonsat' or satellitedata == 'seasonnotsat':

        if satellitedata == 'seasonsat':
            satelliteobserved = ' (satellitobserverad molnfraktion)'
            plt.title(u'Regressionsanalys för alla regioner ' + satelliteobserved + ',' + u'\n medel av månadsmedel (2001-2020)', size=titlesize)

        elif satellitedata == 'seasonnotsat':
            satelliteobserved = ''
            plt.title(u'Regressionsanalys för alla regioner' + satelliteobserved + ',' + u'\n medel av månadsmedel (1979-2020)', size=titlesize)

        #satelliteobserved = ' (satellite observed CF)'
        bestfit, = plt.plot(linexvalues, lineyvalues, linestyle='dashed', color='grey')
        datapointsdjf, = plt.plot(x[::4], y[::4], 'o', color='black', markerfacecolor='black',markersize = datapointsize ,fillstyle = 'full',ms=10)
        datapointsmam, = plt.plot(x[1::4], y[1::4], 'o', color='black', markerfacecolor='black', markersize = datapointsize, fillstyle = 'top',ms=10)
        datapointsjja, = plt.plot(x[2::4], y[2::4], 'o', color='black', markerfacecolor='black',markersize = datapointsize,fillstyle = 'bottom',ms=10)
        datapointsson, = plt.plot(x[3::4], y[3::4], 'o', color='black', markerfacecolor='none', markersize = datapointsize,fillstyle = 'none',ms=10)
        #plt.title('Bestfit for all regions ' + satelliteobserved + '\n mean of monthly means (2001-2020)', size=20)
        plt.title(u'Regressionsanalys för alla regioner' + satelliteobserved + ',' + u'\n medel av säsongsmedel (1979-2020)',size=titlesize)
        plt.ylim(min(y) - 0.1, max(y) + 0.1)
        plt.legend([datapointsdjf,datapointsmam,datapointsjja,datapointsson,bestfit], legendlist, fontsize=legendsize)

        k = 0
        q = 0
        for l in range(20):
            plt.text(x[l]+0.10, y[l], str(regionnames[q]), fontsize=regionsize)
            k += 1
            if k == 4:
                q += 1
                k = 0

    # Creates scatterplot
    plt.xlim(min(x) - 1, max(x) + 1)
    plt.xlabel('LTS [K]', size = xtitlesize)
    plt.ylabel(u'Molnfraktion låga moln', size = ytitlesize)
    plt.xticks(size=xtickssize)
    plt.yticks(size=ytickssize)
    plt.show()

    if satellitedata == 'seasonsat' or satellitedata == 'seasonnotsat':
        pass
    else:
        print(regionlist[i] + " a & b are estimated as: \n a = "
              + str(round(popt[0], 3)) + " " + u"\u00B1" + " " + str(round(pcov[0, 0] ** 0.5, 3)) + "\n"
              + " b = " + str(round(popt[1], 2)) + " " + u"\u00B1" + " " + str(round(pcov[1, 1] ** 0.5, 2)))
        print('Correlation coefficient for ' + regionlist[i] + satelliteobserved + ' is: %.3f' % corr)
        print(" ")


#Calculates mean values for CF and LTS and then scatterplots them
def LTSCF_calculations_and_plot(amount_of_plots_meanmonth,amount_of_plots_allmonths,CFmonthvalues,LTSmonthvalues):
    CFmonthvalues = np.array(CFmonthvalues)
    LTSmonthvalues = np.array(LTSmonthvalues)
    monthnames = ['jan','feb','mar','apr','maj','jun','jul','aug','sep','okt','nov','dec']
    global monthnames

    ###Scatterplot for each region, where mean is taken of all monthly means (12 data points)###

    #Loops all regions and plots scatterplot of mean of monthly means of LTS and CF (not satellite observed)
    s = 0
    for i in range(amount_of_plots_meanmonth):
        #Define scatterplott datapoints
        LTSmonthvaluesnotsat = LTSmonthvalues[i]
        CFmonthvaluesnotsat = CFmonthvalues[i]

        #Call scatterplot function
        scatterplot(LTSmonthvaluesnotsat, CFmonthvaluesnotsat, False)

        #Condition for scatterplot function
        s += 1

    ### Scatterplot for each region, where datapoints are monthly mean measurement of LTS and CF(datapoint amount = 12 * years) ###
    MeanT1000values = []
    MeanT700values = []
    MeanCFvalues = []

    #Append monthly mean temperature values in each region at 1000 and 700 hPa and CF values into lists above.
    for i in range(len(regionlist)):
        MeanT1000values.append(T1000values[i].mean(-1).mean(-1))
        MeanT700values.append(T700values[i].mean(-1).mean(-1))
        MeanCFvalues.append(CFvalues[i].mean(-1).mean(-1))

    #Convert the lists into arrays
    MeanT1000values = np.array(MeanT1000values)
    MeanT700values = np.array(MeanT700values)
    MeanCFvalues = np.array(MeanCFvalues)

    #Calculate LTS for each month
    for i in range(len(MeanT1000values)):
        Pot1000scatter = MeanT1000values - 273.15
        Pot700scatter = (MeanT700values * (float(1000) / 700) ** 0.286) - 273.15
        LTSscatter = Pot700scatter - Pot1000scatter

    #Condition for scatterplot function below
    s = 0

    #Plots scatterplot for each region all datapoints, CF not satellite observed
    for i in range(amount_of_plots_allmonths):
        # Define scatterplot datapoints
        LTSmonthvaluesnotsat = LTSscatter[i]
        CFmonthvaluesnotsat = MeanCFvalues[i]

        #Call scatterplot function
        scatterplot(LTSmonthvaluesnotsat, CFmonthvaluesnotsat,'allmonthsnotsat')
        s += 1
    global s, LTSscatter


#Calculates mean values for CF and LTS and then scatterplots them (where CF is from satellite observations)
def LTSCFsat_calculations_and_plot(amount_of_plots_allmonths_sat,amount_of_plots_meanmonth_sat):
    CFsatallmonths = []
    LTSsatseason = []
    CFsatseason = []
    h = 0
    for i in range(5):

        #Read nc file
        dataplevel= readfilepslevel("Pressurelvl_79-20_" + regionlist[i])
        datasatellite = readfilepslevel("Satellit_200101-202012_" + regionlist[i])

        #Extract variable from nc file. Slice temperaturedata so its starts from 2001.
        T1000sat = dataplevel.variables['t'][264:, 1, :, :]  # 4 axes: time, pressure level, lat,lon
        T700sat = dataplevel.variables['t'][264:, 0, :, :]
        CFsat = datasatellite.variables['cldarea_total_daynight_mon'][:,:,]
        CFsatallmonths.append(CFsat.mean(-1).mean(-1)/float(100))

        #Define empty lists for calculations below
        allaar1000sat = []
        allaar700sat = []
        LTSmonthsat = []
        allaarCFsat = []
        allaarCFmonthsat = []

        #Go through every unique month
        for i in range(0, 12):

            #Create empty list for every month
            uniktarslista1000sat = []
            uniktarslista700sat = []
            uniktarslistaCFsat = []

            # k = month of the year
            k = i

            # Go through every year. len(T1000)/12 = amount of years.
            for j in range(0, len(T1000sat) / 12):

                #Append for every january e.g. and so forth to list above
                uniktarslista1000sat.append(T1000sat[k])
                uniktarslista700sat.append(T700sat[k])
                uniktarslistaCFsat.append(CFsat[k])

                #Go 12 steps forward to get to next years january
                k += 12

            #Append all monthly values to defined lists in the beginning
            allaar1000sat.append(uniktarslista1000sat)
            allaar700sat.append(uniktarslista700sat)
            allaarCFsat.append(uniktarslistaCFsat)

        allaar1000sat = np.array(allaar1000sat)
        allaar700sat = np.array(allaar700sat)
        allaarCFsat = np.array(allaarCFsat)


        #Calculate mean monthly mean values for whole time period
        for l in range(0, 12):
            #Calculate mean temperature at 1000 and 700 hPa
            allaar1000m = allaar1000sat[l].mean(0).mean(0).mean(0)
            allaar700m = allaar700sat[l].mean(0).mean(0).mean(0)

            #Calculate mean CF (satellite observed) for each month
            allaarCFm = allaarCFsat[l].mean(0).mean(0).mean(0)

            # Calculate mean LTS per month
            PotT1000 = allaar1000m - 273.15
            PotT700 = (allaar700m * (float(1000) / 700) ** 0.286) - 273.15
            LTS = PotT700 - PotT1000

            #Append these values to empty defined lists
            allaarCFmonthsat.append(allaarCFm)
            LTSmonthsat.append(LTS)

        #Append all monthly means of all regions to lists to calculate seasonal means
        CFsatseason.append(allaarCFmonthsat)
        LTSsatseason.append(LTSmonthsat)

    # Turn these lists into arrays and divide CF by 100 because its in precentages
    LTSsatmeanmonth = np.array(LTSsatseason)
    CFsatmeanmonth = np.array(CFsatseason) / float(100)
    h = 0
    for i in range(amount_of_plots_meanmonth_sat):
        # Call scatterplot function. h is a
        scatterplot(LTSsatmeanmonth[i], CFsatmeanmonth[i], True)
        h += 1

    CFsatallmonths2 = []
    h = 0
    for i in range(5):

        #Read nc file
        datasatellite = readfilepslevel("Satellit_200101-202012_" + regionlist[i])
        CFsat = datasatellite.variables['cldarea_total_daynight_mon'][:,:,:]
        CFsatallmonths2.append(CFsat.mean(-1).mean(-1)/float(100))

    CFsatallmonths2 = np.array(CFsatallmonths2)

    for i in range(amount_of_plots_allmonths_sat):
        scatterplot(LTSscatter[i][264:],CFsatallmonths2[i], 'allmonthssat')
        h += 1

    global h, CFsatmeanmonth, LTSsatmeanmonth




################ Implementation of functions ########################

######################## 1. Calculate LTS/CF for each region #####################
a = "Pressurelvl_79-20_"
b = "Singlelvl_79-20_"
regionlist = ['Peruvian','Namibian','Californian','Australian','Canarian']

filenamelistpressure = []
filenamelistsingle = []

#Creates list of file names
for i in range(len(regionlist)):
    filenamelistpressure.append(a + regionlist[i])
    filenamelistsingle.append(b + regionlist[i])

CFmonthvalues = []
LTSmonthvalues = []

T1000values = []
T700values = []
CFvalues = []
LTSseasonnotsat = []
CFseasonnotsat = []

#Calculates LTS/CF in each KH-region
for i in range(len(regionlist)):
    print("--------------DATA FOR " + regionlist[i] + "---------------")
    LTSCF(filenamelistpressure[i],regionlist[i], True, True) #Calc annual LTS
    LTSCF(filenamelistpressure[i], regionlist[i], True, False) #Calc monthly, djf.. LTS
    print(" ")
    LTSCF(filenamelistsingle[i], regionlist[i], False, True) #Calc annual CF
    LTSCF(filenamelistsingle[i], regionlist[i], False, False) #Calc monthly, djf.. CF
    print(" ")
    print(" ")
    CFmonthvalues.append(allaarCFmonth)
    LTSmonthvalues.append(LTSmonth)
    T1000values.append(T1000)
    T700values.append(T700)
    CFvalues.append(LCC)
    LTSseasonnotsat.append(LTSnotsat)
    CFseasonnotsat.append(CFnotsat)

################# 2. CF(not satellite observed)-LTS Scatterplot #############
LTSCF_calculations_and_plot(1,1,CFmonthvalues,LTSmonthvalues)

################### 3. Scatterplott satellite observed CF ####################

LTSCFsat_calculations_and_plot(1,1)

################### 5. Scatterplott of seasons all regions (CF not satellite observed) ####################
LTSseasonnotsat = np.array(LTSseasonnotsat)
CFseasonnotsat = np.array(CFseasonnotsat)

LTSseasonnotsat = np.hstack(LTSseasonnotsat)
CFseasonnotsat = np.hstack(CFseasonnotsat)
CFseasonnotsat = np.hstack(CFseasonnotsat)

scatterplot(LTSseasonnotsat,CFseasonnotsat, 'seasonnotsat')

################### 4. Scatterplott of seasons all regions (CF satellite observed) ####################

CFseason = CFsatmeanmonth
LTSseason = LTSsatmeanmonth

allCFseasonvalues = []
allLTSseasonvalues = []

for i in range(5):
    k = -1
    CFseasonvalues  = []
    LTSseasonvalues = []

    # Beraknar CF for djf, mam...
    for l in range(4):
        monthlymeansCF = (1 / float(3)) * (CFseason[i][k] + CFseason[i][k + 1]+ CFseason[i][k + 2])
        CFseasonvalues.append(monthlymeansCF)

        monthlymeansLTS = (1 / float(3)) * (LTSseason[i][k] + LTSseason[i][k + 1]+ LTSseason[i][k + 2])
        LTSseasonvalues.append(monthlymeansLTS)
        k += 3

    allCFseasonvalues.append(CFseasonvalues)
    allLTSseasonvalues.append(LTSseasonvalues)

allCFseasonvalues = np.array(allCFseasonvalues)
allLTSseasonvalues = np.array(allLTSseasonvalues)

allCFseasonvalues = np.hstack(allCFseasonvalues)
allLTSseasonvalues = np.hstack(allLTSseasonvalues)

scatterplot(allLTSseasonvalues,allCFseasonvalues, 'seasonsat')

################# 5. Global distribution of low- and total cloud cover ################################
#globaldistrib("Singlelvl_17-20_Global")

################# 6. Timeseriesplot of LTS ######################################
timeseries(regionlist,0,41)