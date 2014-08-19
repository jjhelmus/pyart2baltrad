"""
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
"""
## Mapping of RAVE to/from Py-ART

## @file
## @author Daniel Michelson (SMHI) and Jonathan Helmus (ARM)
## @date 2014-08-19

import time
import datetime

import numpy as np

try:
    from pyart.config import FileMetadata, get_fillvalue
    from pyart.core.radar import Radar
    from pyart.io.common import make_time_unit_str
except ImportError:
    raise ImportError("Can't find Py-ART. Please install it.")
import netCDF4

try:
    import _rave
    import _raveio
    import _projection
    import _polarvolume
    import _polarscan
    import _polarscanparam
    from Proj import dr, rd
except ImportError:
    raise ImportError("Can't find RAVE. Please install it.")


## Takes a Py-ART radar object and maps it to a RAVE polar scan or volume
# @param radar Py-ART radar object
# @return RaveIO object
def radar2raveio(radar):
    """ Map a Py-ART Radar object to a RAVE polar scan or volume """
    # If the radar object contains a single scan, then we return a _polarscan
    # otherwise a _polarvolume
    if radar.nsweeps > 1:
        obj = _polarvolume.new()
    else:
        obj = _polarscan.new()

    # Set top-level attributes
    if (radar.metadata is not None) and ('source' in radar.metadata):
        obj.source = str(radar.metadata["source"])
    # "nominal" time (the time the radar started acqusition)
    # Py-ART Radar object do not store this, use the start time.
    dt_start = datetime.datetime.strptime(
        radar.time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
    obj.date = dt_start.strftime('%Y%m%d')
    obj.time = dt_start.strftime('%H%M%S')
    obj.height = radar.altitude["data"][0]
    obj.longitude = radar.longitude["data"][0] * dr
    obj.latitude = radar.latitude["data"][0] * dr
    obj.beamwidth = (
        radar.instrument_parameters['radar_beam_width_h']['data'][0] * dr)

    # record data in each scan
    if radar.nsweeps == 1:
        _fillscan(obj, radar)
    else:
        for i in range(radar.nsweeps):
            scan = _polarscan.new()
            scan.longitude, scan.latitude = obj.longitude, obj.latitude
            scan.height = obj.height
            _fillscan(scan, radar, index=i)
            obj.addScan(scan)
    rio = _raveio.new()
    rio.object = obj
    return rio


## Fills a scan object with data and metadata from the radar object for sweep
# at index.
# @param scan PolarScan object
# @param radar radar object
# @param index int sweep index with base=0
def _fillscan(scan, radar, index=0):
    """ Fill a RAVE scan with data and metadata from a radar object. """

    startray = radar.sweep_start_ray_index['data'][index]
    stopray = radar.sweep_end_ray_index['data'][index]
    sweep_times = radar.time['data'][startray:stopray+1]

    # Dataset-specific 'where'
    scan.elangle = radar.elevation["data"][startray] * dr
    scan.rstart = float(radar.range["meters_to_center_of_first_gate"])
    scan.rscale = float(radar.range["meters_between_gates"])
    scan.a1gate = np.argmin(sweep_times) + startray
    # These are not settable in RAVE
    #scan.nrays = stopray - startray + 1
    #scan.nbins = radar.ngates

    # Dataset-specific 'what'
    dt_start = netCDF4.num2date(sweep_times.min(), radar.time['units'])
    scan.startdate = dt_start.strftime('%Y%m%d')
    scan.starttime = dt_start.strftime('%H%M%S')
    dt_end = netCDF4.num2date(sweep_times.max(), radar.time['units'])
    scan.enddate = dt_end.strftime('%Y%m%d')
    scan.endtime = dt_end.strftime('%H%M%S')

    # Dataset-specific 'how'.
    # Such optional attributes have to be named specifically.
    scan.addAttribute("how/startazA",
                      radar.azimuth["data"][startray:stopray+1])

    # Quantity/parameter-specific 'what'
    # Py-ART delagates any scaling and offset of data to the
    # field 'data' dictionary object, only the 'final' values are available
    # for general purpose use.  In additional all bad/missing/undetected
    # data is indicated by possible masking.
    # RAVE has conventions for scaling/offset and missing data vs undetected
    # data. These are not used here.
    for quant in radar.fields.keys():
        param = _polarscanparam.new()
        param.quantity = str(quant)
        param.gain = 1.0    # See above discussion
        param.offset = 0.0
        param.nodata = get_fillvalue()
        param.undetect = get_fillvalue()
        sweep_data = radar.fields[quant]['data'][startray:stopray+1]
        param.setData(np.ma.filled(sweep_data, get_fillvalue()))
        scan.addParameter(param)

    # Unambiguous velocity (Nyquist interval)
    if radar.instrument_parameters is not None:
        inst_params = radar.instrument_parameters
        if 'nyquist_velocity' in inst_params:
            scan.addAttribute(
                'how/NI',
                float(inst_params['nyquist_velocity']['data'][startray]))

    # Site-specific navigation with PROJ.4. to make the object "transformable"
    scan.projection = _projection.new(
        'longlat',
        'Site-specific longlat projection',
        '+proj=latlong +ellps=WGS84 +datum=WGS84')
    return


## Reads an ODIM_H5 file and returns a Py-ART radar object
# @param rio RaveIO object
# @return radar object
def raveio2radar(rio):
    """ Map a RaveIO object to a Py-ART Radar. """

    # create metadata retrieval object
    filemetadata = FileMetadata('odim_h5')

    # determine some key parameters
    hfile = rio.object

    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        nsweeps = 1
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        nsweeps = hfile.getNumberOfScans()
    else:
        raise TypeError(
            "Unsupported object, only SCANs and PVOLs supported.")

    # latitude, longitude and altitude
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')
    latitude['data'] = np.array([hfile.latitude*rd])
    longitude['data'] = np.array([hfile.longitude*rd])
    altitude['data'] = np.array([hfile.height])

    # metadata
    metadata = filemetadata('metadata')
    metadata['source'] = hfile.source
    metadata['original_container'] = 'odim_h5'

    # Datasets, read with each sweep
    datasets = []
    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        datasets.append(unicode("dataset1"))
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        for i in range(rio.object.getNumberOfScans()):
            datasets.append(unicode("dataset%i" % (i+1)))

    # sweep_start_ray_index, sweep_end_ray_index
    # Not to be confused with where/a1gate!
    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        rays_per_sweep = [hfile.nrays-1]
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        rays_per_sweep = [hfile.getScan(i).nrays for i in range(nsweeps)]

    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')
    sweep_start_ray_index['data'] = np.cumsum(
        np.append([0], rays_per_sweep[:-1])).astype('int32')
    sweep_end_ray_index['data'] = np.cumsum(
        rays_per_sweep).astype('int32') - 1

    # sweep_number
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # sweep_mode
    sweep_mode = filemetadata('sweep_mode')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])

    # scan_type
    scan_type = 'ppi'

    # fixed_angle
    fixed_angle = filemetadata('fixed_angle')
    sweep_el = []
    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        sweep_el = [hfile.elangle*rd]
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        sweep_el = [hfile.getScan(i).elangle * rd for i in range(nsweeps)]
    fixed_angle['data'] = np.array(sweep_el, dtype='float32')

    # elevation
    # A better solution is to use the elevation angles for each ray if
    # available in how/startelA, how/stopelA in ODIM_H5 v2.2
    elevation = filemetadata('elevation')
    elevation['data'] = np.repeat(sweep_el, rays_per_sweep)

    # range
    _range = filemetadata('range')
    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        obj = hfile
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        obj = hfile.getScan(0)
        # Check that gate spacing is constant for all scans in the pvol
        # Py-ART Radar object do not support radar data where the gate spacing
        # is not constant for all radials.
        # Data of this type should raises an exception.
        rscales = np.array([hfile.getScan(i).rscale for i in range(nsweeps)])
        nbins = np.array([hfile.getScan(i).nbins for i in range(nsweeps)])
        if np.any(nbins != nbins[0]) or np.any(rscales[0] != rscales):
            raise TypeError(
                "Py-ART cannot handle volumes containing scans with",
                "different (bin) gate spacings or dimensions.")
    # This is a generalization, but we'll live with it.
    _range['data'] = (np.arange(obj.nbins, dtype='float32') * obj.rscale +
                      obj.rstart)
    _range['meters_to_center_of_first_gate'] = obj.rstart
    _range['meters_between_gates'] = obj.rscale

    # azimuth
    # azimuth angle for all rays collected in the volume
    azimuth = filemetadata('azimuth')
    total_rays = np.sum(rays_per_sweep)
    az_data = np.ones((total_rays, ), dtype='float32')

    nbins = obj.nbins
    # loop over the sweeps, store the starting azimuth angles.
    # an average of the startazA and stopazA would probably be a better
    # estimate, but the discontinuity between 0 and 360 would need to be
    # addressed. This is attempted if startazA is available.
    start = 0

    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        if 'how/startazA' in hfile.getAttributeNames():
            sweep_az = obj.getAttribute('how/startazA')
            # + (obj.beamwidth*rd*0.5)
            sweep_az = np.where(np.greater(sweep_az, 360.0),
                                sweep_az-360.0, sweep_az)
            az_data[start:start+obj.nrays] = sweep_az
        else:
            az_data = np.arange(obj.nrays)+(360./obj.nrays/2)

    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        for s in range(nsweeps):
            scan = hfile.getScan(s)
            if 'how/startazA' in scan.getAttributeNames():
                sweep_az = scan.getAttribute('how/startazA')
                #+(scan.beamwidth*rd*0.5)
                sweep_az = np.where(np.greater(sweep_az, 360.0),
                                    sweep_az-360.0, sweep_az)
                az_data[start:start+scan.nrays] = sweep_az
                start += scan.nrays
            else:
                az_data = np.arange(scan.nrays)+(360./scan.nrays/2)

    azimuth['data'] = az_data

    # time
    # time at which each ray was collected.  Need to define
    # starting time and the units of the data in the 'unit' dictionary
    # element.
    # Since startazT and stopazT do not appear to be present in all files
    # and the startepochs and endepochs attributes appear the same for
    # each sweep, just interpolate between these values.
    # XXX This is does not seem correct.
    # Assuming these are UTC times
    _time = filemetadata('time')

    # Remembering that obj is either a single scan or the
    # first scan in the pvol
    attrnames = obj.getAttributeNames()
    if 'how/startepochs' in attrnames and 'how/stopepochs' in attrnames:
        start_epoch = obj.getAttribute('how/startepochs')
        end_epoch = obj.getAttribute('how/stopepochs')
    else:
        start_epoch = time.mktime(datetime.datetime.strptime(
            obj.startdate+obj.starttime, "%Y%m%d%H%M%S").timetuple())
        end_epoch = time.mktime(datetime.datetime.strptime(
            obj.enddate+obj.endtime, "%Y%m%d%H%M%S").timetuple())
    start_time = datetime.datetime.utcfromtimestamp(start_epoch)
    delta_sec = end_epoch - start_epoch
    _time['units'] = make_time_unit_str(start_time)
    _time['data'] = np.linspace(0, delta_sec, total_rays).astype('float32')

    # fields
    # the radar moments or fields are stored in as a dictionary of
    # dictionaries.  The dictionary for each field, a 'field dictionary'
    # should contain any necessary metadata.  The actual data is stored in
    # the 'data' key as a 2D array of size (nrays, ngates) where nrays is the
    # total number of rays in all sweeps of the volume, and ngate is the
    # number of bins or gates in each radial.
    fields = {}
    if rio.objectType in (_rave.Rave_ObjectType_SCAN,
                          _rave.Rave_ObjectType_PVOL):
        h_field_keys = []
        if rio.objectType is _rave.Rave_ObjectType_SCAN:
            field_names = obj.getParameterNames()
        elif rio.objectType is _rave.Rave_ObjectType_PVOL:
            # This tactic assumes that all quantities are available in all
            # scans of the pvol, and that all quantities are ordered the
            # same way in each scan. This is probably dangerous.
            field_names = hfile.getScan(0).getParameterNames()
        for i in range(len(field_names)):
            h_field_keys.append(unicode('data%i' % (i+1)))  # matching order

        # loop over the fields, create the field dictionary
        for field_name, h_field_key in zip(field_names, h_field_keys):
            # XXX still need to set metadata, some default field metadata can
            # likely be provided form the filemetadata object.
            #field_dic = filemetadata(field_name)
            field_dic = {}
            # Assumes the same dtype for each quantity. Potentially dangerous.
            if rio.objectType is _rave.Rave_ObjectType_SCAN:
                dtype = obj.getParameter(field_name).getData().dtype
            elif rio.objectType is _rave.Rave_ObjectType_PVOL:
                dtype = hfile.getScan(0).getParameter(
                    field_name).getData().dtype
            field_dic['data'] = np.zeros((total_rays, nbins), dtype=dtype)
            start = 0
            # loop over the sweeps, copy data into correct location in
            # data array
            if rio.objectType is _rave.Rave_ObjectType_SCAN:
                sweep_data = obj.getParameter(field_name).getData()
                field_dic['data'][start:start + obj.nrays] = sweep_data[:]
            elif rio.objectType is _rave.Rave_ObjectType_PVOL:
                for i in range(hfile.getNumberOfScans()):
                    scan = hfile.getScan(i)
                    sweep_data = scan.getParameter(field_name).getData()
                    param = scan.getParameter(field_name)
                    scaled_sweep_data = (sweep_data[:] * param.gain +
                                         param.offset)
                    field_dic['data'][start:start + scan.nrays] = (
                        scaled_sweep_data[:])
                    start += scan.nrays
            fields[field_name] = field_dic

    # instrument_parameters
    # this is also a dictionary of dictionaries which contains
    # instrument parameter like wavelength, PRT rate, nyquist velocity, etc.
    # A full list of possible parameters can be found in section 5.1 of
    # the CF/Radial document.
    # prt, prt_mode, unambiguous_range, and nyquist_velocity are the
    # parameters which we try to set in Py-ART although a valid Radar object
    # can be created with fewer or more parameters
    beam_width_h = filemetadata.get_metadata('radar_beam_width_h')
    beam_width_h['data'] = np.array([hfile.beamwidth * rd], dtype='float32')

    instrument_parameters = {'radar_beam_width_h': beam_width_h}

    return Radar(
        _time, _range, fields, metadata, scan_type,
        latitude, longitude, altitude,
        sweep_number, sweep_mode, fixed_angle, sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth, elevation,
        instrument_parameters=instrument_parameters)
