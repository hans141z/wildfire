import pandas as pd
from urllib.request import urlopen
import json

def get_ca_wildfires(output_file, start_date = 2013, end_date = 2022):
    """
       Retrieve the California Wildfire Incidences and saves them to a csv 
       
       Arguments: output_file (str) path to output file
       Returns: None
    """
    if (start_date < 2013):
        start_date = 2013
    if (end_date > 2022):
        end_date = 2022
    
    years = range(start_date, end_date,1)


    # scrape the urls
    fireData = []
    for year in years:
        url = "https://www.fire.ca.gov/incidents/"+ str(year) +"/"
        page= urlopen(url)
        html_bytes = page.read()
        html = html_bytes.decode("utf-8")


        start_index = html.find("create_map([") + len("create_map([")
        end_index = html.find('<div id="esri-map">' )-40
        data = html[start_index:end_index]
        for fire in data.split("},"):
            fireData.append(json.loads(fire+"}"))

    data = pd.DataFrame.from_records(fireData)



    # get the cause of the wildfire
    causes = []
    i = 0
    for url in data['Url'].values:
        i+=1;
        page = urlopen(url)
        html_bytes = page.read()
        html = html_bytes.decode("utf-8")

        start = html.find("<th>Cause")
        if (start>0 and html[start_index:]!= "\n"):
            start_index = html[start:].find("<td>") + start + len("<td>")
            end_index = html[start_index:].find("</td>") + start_index
            cause = html[start_index:end_index]
            causes.append(cause)
        else:
            causes.append("Unknown")


    # data formatting
    data['Cause'] = causes
    data = data[data['Cause'].str.len()<=30]
    data['Year'] = data['StartedDate'].str.split("-", expand = True)[0]


    data = data[(data['Latitude']<=90) & (data['Latitude']>=-90)]
    dates = data['StartedDate'].str.split("-", expand = True)
    dates.columns = ['Year','Month', 'Day']
    data = data.merge(dates[['Month', 'Day']], left_index=True, right_index = True)
    data['AcresBurnedDisplay'] = data['AcresBurnedDisplay'].str.replace(",", "")
    data['AcresBurnedDisplay'] = pd.to_numeric(data['AcresBurnedDisplay'])



    data.to_csv(output_file, index = False)