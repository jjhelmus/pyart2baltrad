'''
Copyright (C) 2014- Swedish Meteorological and Hydrological Institute (SMHI)

This file is part of RAVE.

RAVE is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RAVE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with RAVE.  If not, see <http://www.gnu.org/licenses/>.
'''
## Mapping of RAVE to/from Py-ART

## @file
## @author Daniel Michelson (SMHI) and Jonathan Helmus (ARM)
## @date 2014-05-20

import time, datetime
import pyart
import netCDF4

import _raveio
import _projection
import _polarvolume, _polarscan, _polarscanparam
from Proj import dr, rd


## Takes a Py-ART radar object and maps it to a RAVE polar scan or volume
# @param radar Py-ART radar object
# @return RaveIO object
def fromPyART(radar):
    # If the radar object contains a single scan, then we return a _polarscan
    # otherwise a _polarvolume
    if radar.nsweeps > 1:
        obj = _polarvolume.new()
    else:
        obj = _polarscan.new()

    # Set top-level attributes
    if (radar.metadata is not None) and  ('source' in radar.metadata):
        obj.source = str(radar.metadata["source"])
    dt = datetime.datetime.strptime(radar.time['units'][14:],
                                    '%Y-%m-%dT%H:%M:%SZ')
    # "nominal" time (the time the radar started acqusition)
    obj.date = dt.strftime('%Y%m%d')    # YYYYmmdd
    obj.time =  dt.strftime('%H%M%S')   # HHMMSS

    obj.height = radar.altitude["data"][0]
    obj.longitude = radar.longitude["data"][0] * dr
    obj.latitude = radar.latitude["data"][0] * dr
    obj.beamwidth = (
        radar.instrument_parameters['radar_beam_width_h']['data'][0] * dr)

    if radar.nsweeps == 1:
        fillScan(obj, radar)
    else:
        for i in range(radar.nsweeps):
            scan = _polarscan.new()
            scan.longitude, scan.latitude = obj.longitude, obj.latitude
            scan.height = obj.height
            fillScan(scan, radar, index=i)
            obj.addScan(scan)

    rio = _raveio.new()
    rio.object = obj
    return rio


## Fills a scan object with data and metadata from the radar object for sweep
# at index.
# @param scan PolarScan object
# @param radar radar object
# @param index int sweep index with base=0
def fillScan(scan, radar, index=0):

    startray = radar.sweep_start_ray_index['data'][index]
    stopray = radar.sweep_end_ray_index['data'][index]

    # Dataset-specific 'where'
    scan.elangle = radar.elevation["data"][startray] * dr
    scan.rstart = float(radar.range["meters_to_center_of_first_gate"])
    scan.rscale = float(radar.range["meters_between_gates"])
    scan.a1gate = int(radar.sweep_start_ray_index["data"][index]) # XXX
    #scan.nrays = stopray - startray + 1
    #scan.nbins = radar.ngates

    # Dataset-specific 'what'
    sweep_times = radar.time['data'][startray:stopray+1]
    dt_start = netCDF4.num2date(sweep_times.min(), radar.time['units'])
    dt_end = netCDF4.num2date(sweep_times.max(), radar.time['units'])
    scan.startdate = dt_start.strftime('%Y%m%d')
    scan.starttime = dt_start.strftime('%H%M%S')
    scan.enddate =  dt_end.strftime('%Y%m%d')
    scan.endtime =  dt_end.strftime('%H%M%S')

    # Dataset-specific 'how'. Such optional attributes have to be named specifically.
    scan.addAttribute("how/startazA", radar.azimuth["data"][startray:stopray+1])

    # Quantity/parameter-specific 'what'
    # Note that coefficients weren't addressed in the mapping from RAVE to Py-ART, which they should. Conventions apply to uint8 data.
    for quant in radar.fields.keys():
        param = _polarscanparam.new()
        param.quantity = str(quant)
        param.gain =  1.0 #<-- float. If not applicable then 1.0
        param.offset = 0.0 #<-- float. If not applicable then 0.0
        param.nodata = pyart.config.get_fillvalue()
            #<-- float. Unradiated areas. Convention: 255.0 for reflectivity
        param.undetect = pyart.config.get_fillvalue()
            #<-- float. Radiated but no echo. Convention: 0.0 for reflectivity

        sweep_data = radar.fields[quant]['data'][startray:stopray+1]
        param.setData(sweep_data)  # Double-check!

        scan.addParameter(param)

    # instrument parameters
    if radar.instrument_parameters is not None:
        if 'nyquist_velocity' in radar.instrument_parameters:
            ray = radar.sweep_start_ray_index['data'][index]
            nyquist_velocity = (
                radar.instrument_parameters['nyquist_velocity']['data'][ray])
            scan.addAttribute('how/NI', float(nyquist_velocity))

    # Site-specific navigation with PROJ.4. Should make the object "transformable"
    scan.projection = _projection.new(
        'longlat',
        'Site-specific longlat projection',
        '+proj=latlong +ellps=WGS84 +datum=WGS84')



if __name__ == "__main__":
    pass
