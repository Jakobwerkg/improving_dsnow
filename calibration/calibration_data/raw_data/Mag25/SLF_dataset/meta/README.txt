# -----------------------------------------------------------------------------------------
Title: Evaluating methods to estimate the water equivalent of new snow from daily snow depth recordings

Date: November 2024

Authors: jan.magnusson@slf.ch; bertrand.cluzet@slf.ch; louis.queno@slf.ch; mott@slf;
  moritz.oberrauch@slf.ch; giulia.mazzotti@slf.ch; marty@slf.ch; jonas@slf.ch

Contact: jan.magnusson@slf.ch

Summary: Measurement from sites with observations of new snow water equivalent used in the publication:
"Magnusson J., Cluzet B., Quéno L., Mott R., Oberrauch M., Mazzotti G., Marty C., Jonas T., 2025, Evaluating
methods to estimate the water equivalent of new snow from daily snow depth recordings, Cold Regions Science
and Technology, https://doi.org/10.1016/j.coldregions.2025.104435". The dataset covers the period from 1st September 2016 
to 31 August 2022, and contains data from sites in Switzerland. All data files include sufficient metadata to 
be self-explanatory. Please refer to the manuscript for information about the various observations contained in this dataset.

Coordinate system: All data are provided in CH1903/LV03 (EPSG:21781). For transforming these coordinates, 
please see online converter https://www.swisstopo.admin.ch/en/coordinates-conversion-navref 
or code in various languages available at https://github.com/ValentinMinder/Swisstopo-WGS84-LV03.

Time information: All timestamps are given in Central European Time (UTC+01:00). For all data related to new 
snow measurements, the timestamps indicate the end of the observation period (i.e., when the measurements
were performed).

Missing data: No-data is denoted by -9999 in all files.

# -----------------------------------------------------------------------------------------

Files and folder content:

- STATION_LIST.txt: Contains names, coordinates and altitudes of all stations

- HN: Daily observations of new snow snow depth

- HNW: Daily observations for new snow water equivalent

- HS: Daily and biweekly observations of snow depth

- SWE: Biweekly observations of snow water equivalent

The name of the files describes their content. Note that PROFILE denotes HS/SWE measurements performed by manual profiling, 
while STAKE denotes HS measurements performed at a fixed location using a graded stake.
