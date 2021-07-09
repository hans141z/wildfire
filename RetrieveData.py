import requests
from urllib.request import Request, build_opener, HTTPCookieProcessor, urlopen
import base64
import h5py
import pygmt
import numpy as np
import os
from requests_oauthlib import OAuth1
from datetime import datetime
from bs4 import BeautifulSoup
import re
import rasterio as rio
import netCDF4 as nc
from osgeo import gdal, osr




username = "hfa2105@gmail.com"
password = "Siberianhusky04"
cred = '{0}:{1}'.format(username, password)
cred = base64.b64encode(cred.encode('ascii')).decode('ascii')

#base =  "../../../../../Dropbox/TDI/Capstone/Data/"

            
    
    
""" Soil Moisture """
def format_soil_URL(yyyy, mm, dd):
    if len(str(mm))<2:
        mm = "0"+str(mm)
    if len(str(dd))<2:
        dd = "0"+ str(dd)
    url01 = "https://n5eil01u.ecs.nsidc.org/DP4/SMAP/SPL3SMP_E.004/{}.{}.{}/SMAP_L3_SM_P_E_{}{}{}_R17000_001.h5".format(yyyy,mm,dd,yyyy,mm,dd)
    url02 = "https://n5eil01u.ecs.nsidc.org/DP4/SMAP/SPL3SMP_E.004/{}.{}.{}/SMAP_L3_SM_P_E_{}{}{}_R17000_002.h5".format(yyyy,mm,dd,yyyy,mm,dd)
    
    return url01, url02


def cmr_read_soil_in_chunks(file_object, chunk_size=1024 * 1024):
    """Read a file in chunks using a generator. Default chunk size: 1Mb."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data



def download_soil_data(yyyy, mm, dd, base):
    urls = format_soil_URL(yyyy,mm,dd)
    try:
        req = Request(urls[0])
        req.add_header('Authorization', 'Basic {0}'.format(cred))
        opener = build_opener(HTTPCookieProcessor())
        response = opener.open(req)
    except:
        req = Request(urls[1])
        req.add_header('Authorization', 'Basic {0}'.format(cred))
        opener = build_opener(HTTPCookieProcessor())
        response = opener.open(req)
    length = int(response.headers['content-length'])

    count = 0
    chunk_size = min(max(length, 1), 1024 * 1024)
    outputFile = base +"Soil_Moisture/Soil_Moisture_"+mm+"_"+dd +"_"+ yyyy+".txt"
    with open(outputFile, 'wb') as out_file:
        for data in cmr_read_soil_in_chunks(response, chunk_size=chunk_size):
            out_file.write(data)
    print("Download:     {}-{}-{} completed".format(mm,dd,yyyy))
    return outputFile


def retrieve_soil_moisture_data(yyyy,mm,dd, base):
    """
        Download the soil moisture data
    """
    files = os.listdir(base+"Soil_Moisture/")
    fname = "Soil_Moisture_"+mm+"_"+dd +"_"+ yyyy+".txt"
    if fname in files:
        print("{} already exists".format(fname))
        return fname
    else:
        download_soil_data(yyyy,mm,dd, base)
        
def extract_soil_data(yyyy,mm,dd, base):
    file = "Soil_Moisture_"+mm+"_"+dd +"_"+ yyyy+".txt"
    if file not in os.listdir(base+"Soil_Moisture/"):
        print("Dataset does not exist. Download the data.")
        return [],[],[]
    else:
        f = h5py.File(base+"Soil_Moisture/"+file,"r")

        group_id = list(f.keys())[2]   # key= Soil_Moisture_Retrieval_Data_PM

        sm_data = f[group_id]["soil_moisture_pm"][:, :]
        # find the null filler value, replace it with np.nan
        nullValue = f[group_id]["soil_moisture_pm"].attrs['_FillValue']

        sm_data[sm_data == nullValue] = np.nan
        sm_data[sm_data == 1] = np.nan



        lat_data = f[group_id]["latitude_centroid_pm"][:,:]
        lat_data[lat_data == nullValue] = np.nan
        lat_data[sm_data == 1] = np.nan


        lon_data = f[group_id]["longitude_centroid_pm"][:,:]
        lon_data[lon_data == nullValue] = np.nan
        lon_data[sm_data == 1] = np.nan


        return lon_data, lat_data, sm_data
    

    
"""  Precipitation   """
def format_precipitation_URL(yyyy,mm,dd):
    url = "https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDL.06/{0}/{1}/3B-DAY-L.MS.MRG.3IMERG.{0}{1}{2}-S000000-E235959.V06.nc4".format(yyyy,mm,dd)
    
    return url

def download_precipitation_data(yyyy,mm,dd, base):
    url = format_precipitation_URL(yyyy,mm,dd)
    fname = base+"Precipitation_GPM/Precipitation_GPM_{}_{}_{}.nc4".format(mm,dd,yyyy)
    
    
    r = requests.get(url)#, auth = auth)
    
    
    
    if (r.status_code == 200):
        open(fname, "wb").write(r.content) 
        print("Download:    {}-{}-{} completed".format(mm,dd,yyyy))
    else:
        print("Failed to download {}-{}-{}. {}".format(mm,dd,yyyy,r.status_code))
        print(url)
    return fname

def retrieve_precipitation_data(yyyy,mm,dd, base):
    if len(str(mm))<2:
        mm = "0"+ str(mm)
    if len(str(dd))<2:
        dd = "0"+str(dd)
    
    files = os.listdir(base+"Precipitation_GPM/")
    fname = "Precipitation_GPM_{}_{}_{}.nc4".format(mm,dd,yyyy)

    if fname in files:
        print("{} already exists".format(fname))
        return fname
    else:
        download_precipitation_data(yyyy,mm,dd, base)
    
    
def extract_precipitation_data(yyyy,mm,dd, base):
    if len(str(mm))<2:
        mm = "0" + str(mm)
    if len(str(dd))<2:
        dd = "0"+dd
        
        
        
    file = "Precipitation_GPM_"+mm+"_"+dd +"_"+ yyyy+".nc4"
    
    if file not in os.listdir(base+"Precipitation_GPM/"):
        print("Dataset does not exist. Download the data.")
        return [],[],[]
    else:
        file = base+"Precipitation_GPM/Precipitation_GPM_{}_{}_{}.nc4".format(mm,dd,yyyy)
        SFLats, SFLons = [36.8, 38.4], [-123.17,-120.37]

        data = nc.Dataset(file)

        xDim, yDim = 3600, 1800
        lats = np.linspace(-90, 91, yDim)
        lons = np.linspace(-180, 181, xDim)

        latInd = np.argwhere((lats>=SFLats[0]) & (lats<=SFLats[1])).flatten()
        lats= lats[latInd[0]: latInd[-1]]
        sfLatInd = [latInd.min(), latInd.max()]
        lonInd = np.argwhere((lons>=SFLons[0]) & (lons<=SFLons[1])).flatten()
        lons = lons[lonInd[0]:lonInd[-1]]
        sfLonInd = [lonInd.min(), lonInd.max()]

        x, y = np.meshgrid(lons, lats)

        dataset = data['precipitationCal'][0,:,:].transpose()
        dataset=dataset[sfLatInd[0]:sfLatInd[1], sfLonInd[0]:sfLonInd[1]]
        fillValue = -9999.9
        dataset[dataset==fillValue] = np.nan


        lons, lats, data = x.flatten(), y.flatten(), dataset.flatten()
        
        return lons, lats, data

    
    
    
    
    

""" Surface Temperature"""

def format_surface_temperature_URL(yyyy,mm,dd):
    terraFolder = "https://e4ftl01.cr.usgs.gov/MOLT/MOD11A1.006/"
    
    if len(str(mm)) < 2:
        mm = "0"+str(mm)
    if len(str(dd)) < 2:
        dd = "0"+str(dd)

    yrDay = str(datetime(int(yyyy), int(mm), int(dd)).timetuple().tm_yday)
    if (len(yrDay)<2):
        yrDay = "00"+str(yrDay)
    elif (len(yrDay)<3):
        yrDay = "0"+str(yrDay)
    else:
        yrDay = str(yrDay)

    stFolder = terraFolder+"{0}.{1}.{2}/".format(yyyy,mm,dd,yrDay)

    r = requests.get(stFolder)
    soup = BeautifulSoup(r.text, 'html.parser')
    regionFile = "^MOD11A1.A{}{}.h08v05.006".format(yyyy,yrDay)

    regionFile = soup.find_all(text=re.compile(regionFile))[0]


    url = stFolder+regionFile
    return url

def download_surface_temperature_data(yyyy,mm,dd, base):
    url = format_surface_temperature_URL(yyyy,mm,dd)
    fname = base+"Surface_Temperature/Surface_Temperature_{}_{}_{}.hdf".format(mm,dd,yyyy)
    
    
    r = requests.get(url)#, auth = auth)
    
    
    
    if (r.status_code == 200):
        open(fname, "wb").write(r.content) 
        print("Download:    {}-{}-{} completed".format(mm,dd,yyyy))
    else:
        print("Failed to download {}-{}-{}. {}".format(mm,dd,yyyy,r.status_code))
        print(url)
    return fname

def retrieve_surface_temperature_data(yyyy,mm,dd, base):
    files = os.listdir(base+"Surface_Temperature/")
    fname = "Surface_Temperature_{}_{}_{}.hdf".format(mm,dd,yyyy)

    if fname in files:
        print("{} already exists".format(fname))
        return fname
    else:
        download_surface_temperature_data(yyyy,mm,dd, base)
        
def extract_surface_temperature_data(yyyy,mm,dd, base):
    fName = base+"Surface_Temperature/Surface_Temperature_{}_{}_{}.hdf".format(mm,dd,yyyy)

    with rio.open(fName) as dataset:
        for name in dataset.subdatasets:
            if (re.search("LST_Day_1km", name)):
                file = name

    subdata = rio.open(file)
    fillValue = subdata.profile['nodata']

    lstData = np.where(subdata.read(1) == fillValue, np.nan, subdata.read(1))

    SFLats, SFLons = [36.8, 38.4], [-123.17,-120.37]
    
    # from     ../Data/ModisSinusoidalGridCoordinates.txt
    ul_lat, ul_lon = 30, -130.5407 
    lr_lat, lr_lon = 40, -103.9134 
    
    nx, ny = lstData.shape
    la = np.linspace(ul_lat, lr_lat, nx)
    lo = np.linspace(ul_lon, lr_lon, ny)

    lon_mod, lat_mod = np.meshgrid(lo, la)
    lat_mod = lat_mod[::-1]

    #find indices where latitudes are within SF lats
    lat_cells = np.where((lat_mod[:,0]>SFLats[0]) & (lat_mod[:,0]<SFLats[1]))[0]

    #find indices where longitudes are within SF lons
    lon_cells = np.where((lon_mod[0]>SFLons[0]) & (lon_mod[0]<SFLons[1]))[0]

    # get the soil moisture data for those cells
    #sfSurfaceTemp = lstData[lat_cells[0]:lat_cells[-1], lon_cells[0]:lon_cells[-1]]



    lons= lon_mod[lat_cells[0]:lat_cells[-1], lon_cells[0]:lon_cells[-1]].flatten()
    lats = lat_mod[lat_cells[0]:lat_cells[-1], lon_cells[0]:lon_cells[-1]].flatten()
    data = lstData[lat_cells[0]:lat_cells[-1], lon_cells[0]:lon_cells[-1]].flatten()
    #convert data to celsius (scaling factor is .02 https://icess.eri.ucsb.edu/modis/LstUsrGuide/usrguide_mod11.html#sds)
    scalingFactor = .02
    data = data*scalingFactor -273.15
    
    
    return lons, lats, data  

""" Land Cover """

def extract_land_cover_data(file):
    raster = gdal.Open(file)
    
    srcband = raster.GetRasterBand(1)
    data = srcband.ReadAsArray(0,0,raster.RasterXSize, raster.RasterYSize)

    # relationship between raster positions (pxl/coordinate) and georeferenced coordinates
    xmin, xpixel, _, ymax, _, ypixel = raster.GetGeoTransform()
    width, height = raster.RasterXSize, raster.RasterYSize
    xmax = xmin+width*xpixel
    ymin = ymax+height*ypixel

    coords = (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)

    # old spatial reference, new spatial reference, transformer
    src_srs = osr.SpatialReference()
    src_srs.ImportFromWkt(raster.GetProjection())

    tgt_srs = osr.SpatialReference()
    tgt_srs.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(src_srs, tgt_srs)


    # transforming the coordinates from the old spatial reference to the new spatial reference
    trans_coords = []
    for x0, y0 in coords:
        x,y,z = transform.TransformPoint(x0,y0)
        trans_coords.append([x,y])

    ul, ur, ll, lr = trans_coords[0], trans_coords[1],trans_coords[3],trans_coords[2]
    
    
    # create latitude and longitude grid
    start = np.linspace(ul[1], ll[1], raster.RasterYSize)
    end = np.linspace(ur[1], lr[1], raster.RasterYSize)

    lon_grid = []
    for s, e in zip(start, end):
        lon_grid.append(np.linspace(s, e, raster.RasterXSize))


    #latitudes  
    start = np.linspace(ul[0], ll[0], raster.RasterYSize)
    end = np.linspace(ur[0], lr[0], raster.RasterYSize)

    lat_grid = []
    for s, e in zip(start, end):
        lat_grid.append(np.linspace(s, e, raster.RasterXSize))

    lon_grid, lat_grid = np.array(lon_grid), np.array(lat_grid)
    return lon_grid, lat_grid, data



def find_nearest_land_cover_datapoint(lon, lat, dataset):
    delta = (30/111000)/2
    # indices within 30 m of lat bounds
    lat_indices = np.where((dataset[1]>lat-delta) & (dataset[1]< lat+delta))
    lat_indices = [(i, j) for i, j in zip(lat_indices[0], lat_indices[1])]

    #indices within 30m of lon bounds
    lon_indices = np.where((dataset[0]> lon-delta) & (dataset[0]< lon+delta))
    lon_indices = [(i, j) for i, j in zip(lon_indices[0], lon_indices[1])]

    lat_indices = set(lat_indices)
    lons_indices = set(lon_indices)
    intersection= lat_indices.intersection(lon_indices)

    coord = intersection.pop()

    nearest_lon = dataset[0][coord[0],coord[1]]
    nearest_lat = dataset[1][coord[0], coord[1]]
    datapoint = dataset[2][coord[0], coord[1]]
    
    return nearest_lon, nearest_lat, datapoint


def visual_check_nearest_land_cover_datapoint(lon, lat, dataset, title):
    # receive an error unless I import then relad the import in this function
    import pygmt 
    import imp 
    imp.reload(pygmt)
    
    delta = (30/111000)
    nearest_lon, nearest_lat, datapoint = find_nearest_land_cover_datapoint(lon, lat, dataset) 
    n = 100
    SFLats, SFLons = [lat-n*delta, lat+n*delta], [lon-n*delta,lon+n*delta]

    fig = pygmt.Figure()
    fig.coast(region=[SFLons[0], SFLons[1], SFLats[0], SFLats[1]],frame = "afg", projection = "N12c" , borders = "2/1p,black", land = "beige", water = "skyblue")
    pygmt.makecpt(cmap = "oleron", series = [dataset[2].min(), dataset[2].max()])
    fig.plot(x =lon, y = lat,  color= "black", style= "c0.2cm", pen = "black")

    fig.plot(x =nearest_lon, y = nearest_lat,  color= "red", style= "c0.2cm", pen = "black")
    fig.colorbar(frame = 'af+l" {}"'.format(title))
    return fig










""" Water Vapor  """

















""" TEST """
def test_everything(yyyy,mm,dd):
    
    #soil moisture
    print("\t\tSOIL MOISTURE\t\t")
    url = format_soil_URL(yyyy,mm,dd)
    print("\nDownloading the following url: ", url[0])
    
    outputFile =download_soil_data(yyyy,mm,dd)
    print("\nData downloaded here: ", outputFile)
    
    outputFile2 = retrieve_soil_moisture_data(yyyy, mm,dd)
    print("Data retrieved here: ", outputFile2)
    
    
    lon,lat, data = extract_soil_data(yyyy,mm,dd)
    print("\nFirst 5 Lons, lats, and data extracted\n\t{}\n\t{}\n\t{}".format(lon.flatten()[0:5],\
                                                                           lat.flatten()[0:5], \
                                                                           data.flatten()[0:5]) )
   
    #precipitation
    print("\t\tPrecipitation\t\t")
    url = format_precipitation_URL(yyyy,mm,dd)
    print("\nDownloading the following url: ", url)

    outputFile = download_precipitation_data(yyyy, mm,dd)
    print("\nData downloaded here: ", outputFile)

    outputFile2 = retrieve_precipitation_data(yyyy,mm,dd)
    print("Data retrieved here: ", outputFile2)

    
    lon, lat, data = extract_precipitation_data(yyyy, mm, dd)
    print("\nFirst 5 Lons, lats, and data extracted\n\t{}\n\t{}\n\t{}".format(lon.flatten()[0:5],\
                                                                           lat.flatten()[0:5], \
                                                                           data.flatten()[0:5]) )     
    
    
    #surface temperature
    print("\t\tSURFACE TEMPERATURE\t\t")

    url = format_surface_temperature_URL(yyyy, mm, dd)
    print("\nDownloading the following url: ", url)

    outputFile = download_surface_temperature_data(yyyy, mm, dd)
    print("\nData downloaded here: ", outputFile)

    outputFile2 = retrieve_surface_temperature_data(yyyy,mm,dd)
    print("Data retrieved here: ", outputFile2)

    lon, lat, data = extract_surface_temperature_data(yyyy, mm, dd)
    print("\nFirst 5 Lons, lats, and data extracted\n\t{}\n\t{}\n\t{}".format(lon.flatten()[0:5],\
                                                                           lat.flatten()[0:5], \
                                                                           data.flatten()[0:5]) )  
    


