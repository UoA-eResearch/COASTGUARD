#!/usr/bin/env python3

import os
from glob import glob
import warnings
import pandas as pd
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map
warnings.filterwarnings("ignore")
from datetime import datetime
from Toolshed import (
    Download,
    Toolbox,
    VegetationLine,
    Transects,
)
import ee
import geopandas as gpd
import time

start = time.time()

ee.Initialize()
ee.Authenticate()

shorelines = gpd.read_file("../CoastSat/shorelines.geojson")
shorelines = shorelines[shorelines.id.str.startswith("nzd")]
transects = gpd.read_file("../CoastSat/transects_extended.geojson")

print(f"{time.time() - start}: Reference shorelines and transects loaded")

for sitename in tqdm(shorelines.id.unique()):
  shorelines[shorelines.id == sitename].to_file(f"Data/referenceLines/{sitename}.geojson")

  dates = ["1900-01-01", "2030-01-01"]
  sat_list = ["L5", "L7", "L8", "L9"]

  cloud_thresh = 0.5
  wetdry = False
  referenceLineShp = sitename + ".geojson"
  max_dist_ref = 200

  filepath = Toolbox.CreateFileStructure(sitename, sat_list)

  # Return AOI from reference line bounding box and save AOI folium map HTML in sitename directory
  referenceLinePath = os.path.join(filepath, "referenceLines", referenceLineShp)
  referenceLineDF = gpd.read_file(referenceLinePath)
  polygon, point, lonmin, lonmax, latmin, latmax = Toolbox.AOIfromLine(
      referenceLinePath, max_dist_ref, sitename
  )

  # It's recommended to convert the polygon to the smallest rectangle (sides parallel to coordinate axes)
  polygon = Toolbox.smallest_rectangle(polygon)

  if len(dates) > 2:
      daterange = "no"
  else:
      daterange = "yes"
  years = list(
      Toolbox.daterange(
          datetime.strptime(dates[0], "%Y-%m-%d"),
          datetime.strptime(dates[-1], "%Y-%m-%d"),
      )
  )

  inputs = {
      "polygon": polygon,
      "dates": dates,
      "daterange": daterange,
      "sat_list": sat_list,
      "sitename": sitename,
      "filepath": filepath,
  }

  inputs = Download.check_images_available(inputs)

  Sat = Download.RetrieveImages(inputs, SLC=False)
  metadata = Download.CollectMetadata(inputs, Sat)

  LinesPath = "Data/" + sitename + "/lines"

  os.makedirs(LinesPath, exist_ok=True)

  projection_epsg, _ = Toolbox.FindUTM(polygon[0][0][1], polygon[0][0][0])

  settings = {
      # general parameters:
      "cloud_thresh": cloud_thresh,  # threshold on maximum cloud cover
      "output_epsg": projection_epsg,  # epsg code of spatial reference system desired for the output
      "wetdry": wetdry,  # extract wet-dry boundary as well as veg
      # quality control:
      "check_detection": False,  # if True, shows each shoreline detection to the user for validation
      "adjust_detection": False,  # if True, allows user to adjust the postion of each shoreline by changing the threhold
      "save_figure": True,  # if True, saves a figure showing the mapped shoreline for each image
      # [ONLY FOR ADVANCED USERS] shoreline detection parameters:
      "min_beach_area": 200,  # minimum area (in metres^2) for an object to be labelled as a beach
      "buffer_size": 250,  # radius (in metres) for buffer around sandy pixels considered in the shoreline detection
      "min_length_sl": 500,  # minimum length (in metres) of shoreline perimeter to be valid
      "cloud_mask_issue": False,  # switch this parameter to True if sand pixels are masked (in black) on many images
      # add the inputs defined previously
      "inputs": inputs,
      "projection_epsg": projection_epsg,
      "year_list": years,
  }

  referenceLine, ref_epsg = Toolbox.ProcessRefline(referenceLinePath, settings)

  settings["reference_shoreline"] = referenceLine
  settings["ref_epsg"] = ref_epsg
  settings["max_dist_ref"] = max_dist_ref
  settings["reference_coreg_im"] = None

  output, output_latlon, output_proj = VegetationLine.extract_veglines(
      metadata, settings, polygon, dates
  )
  print(output)

  # Save output veglines
  Toolbox.SaveConvShapefiles(output, LinesPath, sitename, settings["output_epsg"])
  VegBasePath = "Data/" + sitename + "/lines"
  VeglineShp = glob(LinesPath + "/*veglines.shp")
  VeglineGDF = gpd.read_file(VeglineShp[0])
  VeglineGDF = VeglineGDF.sort_values(
      by="dates"
  )  # sort GDF by dates to ensure transect intersects occur in chronological order
  VeglineGDF = VeglineGDF.reset_index(drop=True)  # reset GDF index after date sorting

  TransectGDF = transects[transects.site_id == sitename]
  TransectGDF["reflinepnt"] = transects.geometry.centroid
  TransectGDF.rename(columns={"id": "TransectID"}, inplace=True)
  TransectGDF = TransectGDF[["site_id", "TransectID", "geometry", "reflinepnt"]]

  TransectInterGDF = Transects.GetIntersections(LinesPath, TransectGDF, VeglineGDF)
  dt_lookup = {}
  for _, row in TransectInterGDF.iterrows():
      for i in range(len(row.dates)):
          dt = row.dates[i] + " " + row.times[i]
          if dt not in dt_lookup:
              dt_lookup[dt] = {"dates": dt, "satname": row.satname[i]}
          dt_lookup[dt][row.TransectID] = row.distances[i]

  df = pd.DataFrame(dt_lookup.values())
  fn = os.path.join(
      settings["inputs"]["filepath"],
      settings["inputs"]["sitename"],
      "transect_time_series.csv",
  )
  df.to_csv(fn, index=False, float_format="%.2f")
