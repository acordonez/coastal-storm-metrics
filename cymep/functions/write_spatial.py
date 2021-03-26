import json
import numpy as np
import netCDF4 as nc
from datetime import datetime
import os

output_json = "output.json"

def create_dims(dims_dict, **kwargs):
    """Build dimensions dictionary for json from component dictionaries.
    Parameters:
        dims_dict (dictionary): dictionary to store dimensions in
        **kwargs (name/value pairs): sub-dictionaries for each dimension
    """
    for key, value in kwargs.items():
        dims_dict.update({key: value})
    return dims_dict

def set_metric_definitions(metric_list, desc):
    """Populates metrics definitions using template."""
    metric_dict = {}
    for metric in metric_list:
        pre,suf = metric.split('_')
        metric_desc = desc["prefix"][pre] + desc["suffix"][suf]
        metric_dict[metric] = metric_desc
    return metric_dict

def get_dimensions(json_dict, json_structure):
    """Populate dimensions details from results dictionary.
    Parameters:
        json_dict (dictionary): metrics to pull dimension details from
        json_structure (list): ordered list of dimension names
    """
    keylist = {}
    level = 0
    while level < len(json_structure):
        if isinstance(json_dict, dict):
            first_key = list(json_dict.items())[0][0]
            if first_key == "attributes":
                first_key = list(json_dict.items())[1][0]
            dim = json_structure[level]
            if dim == "statistic":
                keys = [key for key in json_dict]
                keylist[dim] = {"indices": keys}
            else:
                keys = {key: {} for key in json_dict if key != "attributes"}
                keylist[dim] = keys
            json_dict = json_dict[first_key]
        level += 1
    return keylist

def get_env():
    """Return versions of dependencies."""
    import pandas
    import scipy
    import netCDF4
    import subprocess

    versions = {}
    versions['netCDF4'] = netCDF4.__version__
    # numpy is already imported
    versions['numpy'] = np.__version__
    versions['pandas'] = pandas.__version__
    versions['scipy'] = scipy.__version__
    ncl = subprocess.check_output(['ncl', '-V']).decode("utf-8").rstrip()
    versions['ncl'] = ncl
    return versions

def create_output_json(wkdir,test_path,obs_path):

  log_path = wkdir + "/CyMeP.log.txt"
  out_json = {'index': 'index',
              'provenance': {},
              'plots': {},
              'html': "index.html",
              'metrics': {}}
  out_json["provenance"] = {"environment": get_env(),
              "modeldata": test_path,
              "obsdata": obs_path,
              "log": log_path}
  out_json["data"] = {}
  out_json["metrics"] = {}
  out_json["plots"] = {}
  out_json["html"] = {}

  outfilepath = os.path.join(wkdir,output_json)
  if os.path.exists(outfilepath):
    os.remove(outfilepath)

  with open(outfilepath,"w") as outfilename:
    json.dump(out_json, outfilename, indent=2)

def write_spatial_netcdf(spatialdict,permondict,peryrdict,taydict,modelsin,nyears,nmonths,latin,lonin,globaldict,wkdir,cmec=None):

  # Convert modelsin from pandas to list
  modelsin=modelsin.tolist()

  # Set up dimensions
  nmodels=len(modelsin)
  nlats=latin.size
  nlons=lonin.size
  nchar=16

  netcdfdir=wkdir + "/netcdf-files/"
  os.makedirs(os.path.dirname(netcdfdir), exist_ok=True)
  netcdfile=netcdfdir+"/netcdf_"+globaldict['strbasin']+"_"+os.path.splitext(globaldict['csvfilename'])[0]

  # open a netCDF file to write
  ncout = nc.Dataset(netcdfile+".nc", 'w', format='NETCDF4')

  # define axis size
  ncout.createDimension('model', nmodels)  # unlimited
  ncout.createDimension('lat', nlats)
  ncout.createDimension('lon', nlons)
  ncout.createDimension('characters', nchar)
  ncout.createDimension('months', nmonths)
  ncout.createDimension('years', nyears)

  # create latitude axis
  lat = ncout.createVariable('lat', 'f', ('lat'))
  lat.standard_name = 'latitude'
  lat.long_name = 'latitude'
  lat.units = 'degrees_north'
  lat.axis = 'Y'

  # create longitude axis
  lon = ncout.createVariable('lon', 'f', ('lon'))
  lon.standard_name = 'longitude'
  lon.long_name = 'longitude'
  lon.units = 'degrees_east'
  lon.axis = 'X'

  # Write lon + lat
  lon[:] = lonin[:]
  lat[:] = latin[:]

  # create variable arrays
  # Do spatial variables
  for ii in spatialdict:
    vout = ncout.createVariable(ii, 'f', ('model', 'lat', 'lon'), fill_value=1e+20)
   # vout.long_name = 'density'
   # vout.units = '1/year'
    vout[:] = np.ma.masked_invalid(spatialdict[ii][:,:,:])

  # create variable array
  for ii in permondict:
    vout = ncout.createVariable(ii, 'f', ('model', 'months'), fill_value=1e+20)
    vout[:] = np.ma.masked_invalid(permondict[ii][:,:])

  # create variable array
  for ii in peryrdict:
    vout = ncout.createVariable(ii, 'f', ('model', 'years'), fill_value=1e+20)
    vout[:] = np.ma.masked_invalid(peryrdict[ii][:,:])

  # create variable array
  for ii in taydict:
    vout = ncout.createVariable(ii, 'f', ('model'), fill_value=1e+20)
    vout[:] = np.ma.masked_invalid(taydict[ii][:])

  # Write model names to char
  model_names = ncout.createVariable('model_names', 'c', ('model', 'characters'))
  model_names[:] = nc.stringtochar(np.array(modelsin).astype('S16'))

  #today = datetime.today()
  ncout.description = "Coastal metrics processed data"
  ncout.history = "Created " + datetime.today().strftime('%Y-%m-%d-%H:%M:%S')
  for ii in globaldict:
    ncout.setncattr(ii, str(globaldict[ii]))

  # close files
  ncout.close()

  if cmec is not None:
    # Update metadata in output.json
    desc = {
      "region": globaldict["strbasin"],
      "filename": netcdfile,
      "longname": globaldict["strbasin"] + " netcdf output",
      "description": "Coastal metrics processed data"
    }
    with open(os.path.join(wkdir,output_json), "r") as outfilename:
      tmp = json.load(outfilename)
    tmp["data"].update({"netcdf": desc})
    with open(os.path.join(wkdir,output_json), "w") as outfilename:
      json.dump(tmp, outfilename, indent=2)

    # Dump metrics in json format
    write_spatial_jsons(permondict,peryrdict,taydict,modelsin,nyears,globaldict,wkdir)

def write_spatial_jsons(permondict,peryrdict,taydict,modelsin,nyears,globaldict,wkdir):
  # Output metrics from netCDF into JSONs

  with open(os.path.join(wkdir,output_json), "r") as outfilename:
    outjson = json.load(outfilename)

  base_file_name=wkdir+"/json/netcdf_"+globaldict['strbasin']+"_"+os.path.splitext(globaldict['csvfilename'])[0]

  # Load descriptions for all the different file names
  with open("functions/output_templates/output_desc.json") as description_file:
    data_desc = json.load(description_file)

  permon_json_name = base_file_name + "_month.json"
  peryr_json_name = base_file_name + "_year.json"
  tay_json_name = base_file_name + "_taylor.json"
  model_dict = dict.fromkeys(modelsin, {})
  model_count = len(modelsin)

  # Create a descriptive dictionary of these new json files
  # for output.json since these are not described in output_desc.json
  netcdf_desc = {peryr_json_name: {
                                    "filename": "",
                                    "longname": "yearly cyclone metrics",
                                    "description": "yearly cyclone metrics converted from netcdf"
               },
               permon_json_name: {
                                    "filename": "",
                                    "longname": "monthly cyclone metrics",
                                    "description": "monthly cyclone metrics converted from netcdf"
               },
               tay_json_name: {
                                    "filename": "",
                                    "longname": "Taylor cyclone metrics",
                                    "description": "Taylor cyclone metrics converted from netcdf"
               }
  }

  try:
    os.mkdir(wkdir+"/json")
  except FileExistsError:
    pass

  #-----Monthly-----
  # Set up metric dimensions
  json_month = {"DIMENSIONS": {
                "json_structure": ["model", "metric", "month"], "dimensions": {}}}
  metric_list = [*permondict]
  met_dict = set_metric_definitions(metric_list, data_desc["statistics"])
  month_dict = {"0": "January", "1": "February", "2": "March", "3": "April",
                "4": "May", "5": "June", "6": "July", "7": "August",
                "8": "September", "9": "October", "10": "November",
                "11": "December"}
  # Populate RESULTS
  results_json = {"RESULTS": {}}
  for model_num,model_name in enumerate(modelsin):
    results_json["RESULTS"][model_name] = {}
    for metric in metric_list:
      results_json["RESULTS"][model_name].update({metric: {}})
      for time,value in enumerate(permondict[metric][model_num]):
        if np.isnan(value):
          value = None
        results_json["RESULTS"][model_name][metric].update({time: value})
  # Create other json fields
  json_month.update(results_json)
  json_month["DIMENSIONS"]["dimensions"].update({"model":model_dict})
  json_month["DIMENSIONS"]["dimensions"].update({"metric": met_dict})
  json_month["DIMENSIONS"]["dimensions"].update({"month": month_dict})
  # Write json
  with open(permon_json_name, "w") as mfile:
      json.dump(json_month, mfile, indent=2)

  # Update entry in metadata json
  outjson["metrics"][os.path.basename(permon_json_name)] = {
    "longname": "Monthly cyclone metrics",
    "filename": permon_json_name,
    "description": "monthly cyclone metrics converted from netcdf"
  }

  #-----Yearly-----
  json_year = {"DIMENSIONS": {
              "json_structure": ["model", "metric", "year"], "dimensions": {}}}
  metric_list = [*peryrdict]
  met_dict = set_metric_definitions(metric_list, data_desc["statistics"])
  year_list = [y for y in range(int(globaldict["styr"]), int(globaldict["enyr"])+1)]
  # Populate RESULTS
  results_json = {"RESULTS": {}}
  for model_num,model_name in enumerate(modelsin):
    results_json["RESULTS"][model_name] = {}
    for metric in metric_list:
      results_json["RESULTS"][model_name].update({metric: {}})
      for time,value in enumerate(peryrdict[metric][model_num]):
        if np.isnan(value):
          value = None
        time = str(time + int(globaldict["styr"]))
        results_json["RESULTS"][model_name][metric].update({time: value})
  # Create other json fields
  json_year.update(results_json)
  year_dict = dict.fromkeys(year_list, {})
  json_year.update(results_json)
  json_year["DIMENSIONS"]["dimensions"].update({"model":model_dict})
  json_year["DIMENSIONS"]["dimensions"].update({"metric": met_dict})
  json_year["DIMENSIONS"]["dimensions"].update({"year": year_dict})
  # Write json
  with open(peryr_json_name, "w") as yfile:
      json.dump(json_year, yfile, indent=2)

  # Update entry in metadata json
  outjson["metrics"][os.path.basename(peryr_json_name)] = {
    "longname": "Yearly cyclone metrics",
    "filename": peryr_json_name,
    "description": "Yearly cyclone metrics converted from netcdf"
  }

  #-----Taylor-----
  json_taylor = {"DIMENSIONS": {
                "json_structure": ["model", "metric"], "dimensions": {}}}
  metric_list = [*taydict]
  met_dict = set_metric_definitions(metric_list, data_desc["statistics"])
  # Populate Results
  results_json = {"RESULTS": {}}
  for model_num,model_name in enumerate(modelsin):
    results_json["RESULTS"][model_name] = {}
    for metric in metric_list:
      metric_dict = taydict[metric][model_num]
      if np.isnan(metric_dict):
        metric_dict = None
      results_json["RESULTS"][model_name].update({metric: metric_dict})
  # Update other fields
  json_taylor.update(results_json)
  json_taylor["DIMENSIONS"]["dimensions"].update({"model":model_dict})
  json_taylor["DIMENSIONS"]["dimensions"].update({"metric": met_dict})
  # write json
  with open(tay_json_name, "w") as tfile:
    json.dump(json_taylor, tfile, indent=2)

  # Update entry in metadata json
  outjson["metrics"][os.path.basename(tay_json_name)] = {
    "longname": "Taylor cyclone metrics",
    "filename": tay_json_name,
    "description": "Taylor cyclone metrics converted from netcdf"
  }

  # Return metadata info
  with open(os.path.join(wkdir,output_json), "w") as outfilename:
    json.dump(outjson, outfilename, indent=2)

def write_single_json(vardict,modelsin,wkdir,jsonname,statistics,cmec):
  # CSV files
  # This section converts the metrics stored in
  # csv files to the CMEC json format.
  json_structure = ["model", "metric"]
  cmec_json = {"DIMENSIONS": {"json_structure": json_structure}, "RESULTS": {}}
  # convert all the csv files in csv-files/ to json
  if isinstance(modelsin,str):
    modelsin = [modelsin]
  results = dict.fromkeys(modelsin,{})
  if len(modelsin) > 1:
    for ii in vardict:
      for model_num,model in enumerate(modelsin):
        results[model][ii] = vardict[ii][model_num]
      if np.isnan(results[model][ii]):
        results[model][ii] = None
  else:
    for ii in vardict:
      results[modelsin[0]][ii] = vardict[ii]
      if np.isnan(results[modelsin[0]][ii]):
        results[modelsin[0]][ii] = None
  dimensions = get_dimensions(results.copy(), json_structure)
  metric_list = [key for key in dimensions["metric"]]
  dimensions["metric"] = set_metric_definitions(metric_list, statistics)
  cmec_json["DIMENSIONS"]["dimensions"] = dimensions
  cmec_json["RESULTS"] = results

  jsondir = os.path.join(wkdir,"json")
  os.makedirs(jsondir, exist_ok=True)
  cmec_file_name = os.path.join(jsondir,jsonname + ".json")
  with open(cmec_file_name, "w") as cmecfile:
      json.dump(cmec_json, cmecfile, indent=2)

  desc = {cmec[0] + " json": {
    "longname": cmec[1],
    "description": cmec[2],
    "filename": os.path.relpath(cmec_file_name, start=wkdir)
  }}
  with open(os.path.join(wkdir,output_json),"r") as outfilename:
    tmp = json.load(outfilename)
  tmp["metrics"].update(desc)
  with open(os.path.join(wkdir,output_json),"w") as outfilename:
    json.dump(tmp, outfilename, indent=2)

def write_dict_csv(vardict,modelsin):
  # create variable array
  csvdir="./csv-files/"
  os.makedirs(os.path.dirname(csvdir), exist_ok=True)
  for ii in vardict:
    csvfilename = csvdir+"/"+str(ii)+".csv"
    if vardict[ii].shape == modelsin.shape:
      tmp = np.concatenate((np.expand_dims(modelsin, axis=1),np.expand_dims(vardict[ii], axis=1)), axis=1)
    else:
      tmp = np.concatenate((np.expand_dims(modelsin, axis=1), vardict[ii]), axis=1)
    np.savetxt(csvfilename, tmp, delimiter=",", fmt="%s")

def write_single_csv(vardict,modelsin,wkdir,csvname,statistics,cmec=None):
  # create variable array
  csvdir = os.path.join(wkdir,'csv-files')
  os.makedirs(os.path.dirname(csvdir), exist_ok=True)
  csvfilename = csvdir+"/"+csvname+".csv"

  # If a single line CSV with one model
  if np.isscalar(modelsin):
    tmp = np.empty((1,len(vardict)))
    headerstr="Model"
    iterix = 0
    for ii in vardict:
      headerstr=headerstr+","+ii
      tmp[0,iterix]=vardict[ii]
      iterix += 1

    # Create a dummy numpy string array of "labels" with the control name to append as column #1
    labels = np.empty((1,1),dtype="<U10")
    labels[:] = modelsin
    # Stack labels and numpy dict arrays horizontally as non-header data
    tmp = np.hstack((labels, tmp))

  # Else, the more common outcome; 2-D arrays
  else:
    # Concat models to first axis
    firstdict=list(vardict.keys())[0]
    headerstr="Model,"+firstdict

    if vardict[firstdict].shape == modelsin.shape:
      tmp = np.concatenate((np.expand_dims(modelsin, axis=1),np.expand_dims(vardict[firstdict], axis=1)), axis=1)
    else:
      tmp = np.concatenate((np.expand_dims(modelsin, axis=1), vardict[firstdict]), axis=1)

    for ii in vardict:
      if ii != firstdict:
        tmp = np.concatenate((tmp, np.expand_dims(vardict[ii], axis=1)), axis=1)
        headerstr=headerstr+","+ii

  # Write header + data array
  np.savetxt(csvfilename, tmp, delimiter=",", fmt="%s", header=headerstr, comments="")

  if cmec is not None:
    # write to output json
    with open(os.path.join(wkdir,output_json), "r") as outfilename:
      outjson = json.load(outfilename)
    key = cmec[0] + " csv"
    desc = {key: {
      "longname": cmec[1],
      "description": cmec[2],
      "filename": os.path.relpath(csvfilename, start=wkdir)
    }}
    outjson["metrics"].update(desc)
    with open(os.path.join(wkdir,output_json), "w") as outfilename:
      json.dump(outjson, outfilename, indent=2)

    write_single_json(vardict,modelsin,wkdir,csvname,statistics,cmec)
