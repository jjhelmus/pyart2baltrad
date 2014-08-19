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
## @date 2014-04-22

import sys
import numpy as np
import time, datetime

import _rave
from Proj import dr, rd

try:
    from pyart.config import FileMetadata, get_fillvalue
    from pyart.core.radar import Radar
    from pyart.io.common import make_time_unit_str
    from pyart.io.common import radar_coords_to_cart
except ImportError:
    raise ImportError, "Can't find Py-ART. Please install it."


## Reads an ODIM_H5 file and returns a Py-ART radar object
# @param rio RaveIO object
# @return radar object
def toPyART(rio):
    # create metadata retrieval object
    # this object is used to create generic metadata for the various
    # parameters passed to the Radar object.  Using this object allows the
    # user to adjust the default values using the Py-ART config file.
    filemetadata = FileMetadata('odim_h5')

    # determine some key parameters
    hfile = rio.object

    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        nsweeps = 1
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        nsweeps = hfile.getNumberOfScans()
    else:
        raise TypeError, "rave_pyart only supports conversion of SCANs and PVOLs. Exiting..."

    # The general procedure for each parameter is to create a dictionary
    # with default values from the filemetadata object. Then updating any
    # keys with data from the file, as well as recording the actual parameter
    # data in the 'data' key.

    # latitude, longitude and altitude
    # latitude and longitude are measured in degrees (north and east).
    # altitude in meters above mean sea level
    latitude = filemetadata('latitude')
    longitude = filemetadata('longitude')
    altitude = filemetadata('altitude')

    latitude['data'] = np.array([hfile.latitude*rd])
    longitude['data'] = np.array([hfile.longitude*rd])
    altitude['data'] = np.array([hfile.height])

    # metadata
    # this dictionary is used to store any metadata not recorded
    # elsewhere, typically this is quite sparse.
    # The default keys and values are as follow:
    # 'Conventions': 'CF/Radial instrument_parameters',
    # 'comment': '',
    # 'history': '',
    # 'institution': '',
    # 'instrument_name': '',
    # 'references': '',
    # 'source': '',
    # 'title': '',
    # 'version': '1.3'}
    # See section 4.1 of the CF/Radial format for a description of these
    # attributes.  The 'what/source attribute can probably be used to
    # fill in a number of these.
    metadata = filemetadata('metadata')
    metadata['source'] = hfile.source
    metadata['original_container'] = 'odim_h5'

#    # Datasets, read with each sweep
    datasets = []
    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        datasets.append(unicode("dataset1"))
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        for i in range(rio.object.getNumberOfScans()):
            datasets.append(unicode("dataset%i" % (i+1)))

    # sweep_start_ray_index, sweep_end_ray_index
    # These two dictionaries contain the indices of the first and last ray
    # in each sweep (0-based indexing). Not to be confused with where/a1gate!
    sweep_start_ray_index = filemetadata('sweep_start_ray_index')
    sweep_end_ray_index = filemetadata('sweep_end_ray_index')

    rays_per_sweep, ssri, seri = [], [], []
    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        rays_per_sweep.append(hfile.nrays)
        ssri.append(0)
        seri.append(hfile.nrays-1)
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        for s in range(nsweeps):
            scan = hfile.getScan(s)
            rays_per_sweep.append(scan.nrays)
            if len(ssri) == 0:
                ssri.append(0)
                seri.append(scan.nrays-1)
            else:
                ssri.append(ssri[-1] + scan.nrays)
                seri.append(seri[-1] + scan.nrays)

    sweep_start_ray_index['data'] = np.array(ssri, dtype='int32')
    sweep_end_ray_index['data'] = np.array(seri, dtype='int32')

    # sweep_number
    # sweep number in the volume, 0-based
    sweep_number = filemetadata('sweep_number')
    sweep_number['data'] = np.arange(nsweeps, dtype='int32')

    # sweep_mode
    # The scan type for each sweep, common options are:
    # 'rhi','azimuth_surveillance', 'vertical_pointing'
    # Are all ODIM_H5 files PPI volumes? How do you check for this?
    sweep_mode = filemetadata('sweep_mode')
    sweep_mode['data'] = np.array(nsweeps * ['azimuth_surveillance'])

    # scan_type
    # 'ppi', 'rhi', or 'vpt'.
    # Assuming all ODIM_H5 files are ppi scans
    scan_type = 'ppi'

    # fixed_angle
    # this is the elevation or azimuth angle that is fixed for each sweep.
    # In this case the elevation angle.
    fixed_angle = filemetadata('fixed_angle')
    sweep_el = []
    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        sweep_el.append(hfile.elangle*rd)
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        for s in range(nsweeps):
            scan = hfile.getScan(s)
            sweep_el.append(scan.elangle*rd)
    fixed_angle['data'] = np.array(sweep_el, dtype='float32')

    # elevation
    # elevation for each ray in the volume. Since these are PPI scans
    # the elevation angle for each sweep is repeated the number of rays
    # contained in that sweep.
    elevation = filemetadata('elevation')
    elevation['data'] = np.repeat(sweep_el, rays_per_sweep)
    # better solution is to use the elevation angles per ray if available:
    # how/startelA, how/stopelA in ODIM_H5 v2.2

    # range
    # range contains the distances in meters to the center of each range
    # bin.  The 'meters_to_center_of_first_gate' and 'meters_between_gates'
    # attribute should also be set accordingly.  If the gate spacing is not
    # constant, remove the 'meters_beween_gates' key and change
    # 'spacing_is_constant' to 'false'.
    _range = filemetadata('range')

    # here we assume that the gate layout in the first sweep is the same as all
    # the sweeps ("anglesync").  A check should be added to verify this
    # assumption.  The Radar object cannot work with radar data where the
    # gate spacing is not constant for all radials.  Data of this type
    # should raise an exception.
    if rio.objectType is _rave.Rave_ObjectType_SCAN: obj = hfile
    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        # Check that gate spacing is constant for all scans in the pvol
        nbins, rscales = [], []
        for i in range(hfile.getNumberOfScans()):
            scan = hfile.getScan(i)
            rscales.append(scan.rscale)
            nbins.append(scan.nbins)
        if all(x == rscales[0] for x in rscales) and all(b == nbins[0] for b in nbins):
            obj = hfile.getScan(0)
        else:
            raise TypeError, "Py-ART cannot handle volumes containing scans with different (bin) gate spacings or dimensions."
    first_gate = obj.rstart
    gate_spacing = obj.rscale
    nbins = obj.nbins
    # This is a generalization, but we'll live with it.
    _range['data'] = (np.arange(nbins, dtype='float32') * gate_spacing +
                      first_gate)
    _range['meters_to_center_of_first_gate'] = first_gate
    _range['meters_between_gates'] = gate_spacing

    # azimuth
    # azimuth angle for all rays collected in the volume
    azimuth = filemetadata('azimuth')
    total_rays = np.sum(rays_per_sweep)
    az_data = np.ones((total_rays, ), dtype='float32')

    # loop over the sweeps, store the starting azimuth angles.
    # an average of the startazA and stopazA would probably be a better
    # estimate, but the discontinuity between 0 and 360 would need to be
    # addressed. This is attempted if startazA is available.
    start = 0

    if rio.objectType is _rave.Rave_ObjectType_SCAN:
        if 'how/startazA' in hfile.getAttributeNames():
            sweep_az = obj.getAttribute('how/startazA')# + (obj.beamwidth*rd*0.5)
            sweep_az = np.where(np.greater(sweep_az, 360.0),
                                sweep_az-360.0, sweep_az)
            az_data[start:start+obj.nrays] = sweep_az
        else:
            az_data = np.arange(obj.nrays)+(360./obj.nrays/2)

    elif rio.objectType is _rave.Rave_ObjectType_PVOL:
        for s in range(nsweeps):
            scan = hfile.getScan(s)
            if 'how/startazA' in scan.getAttributeNames():
                sweep_az=scan.getAttribute('how/startazA')#+(scan.beamwidth*rd*0.5)
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

    # Remembering that obj is either a single scan or the first scan in the pvol
    attrnames = obj.getAttributeNames()
    if 'how/startepochs' in attrnames and 'how/stopepochs' in attrnames:
        start_epoch = obj.getAttribute('how/startepochs')
        end_epoch = obj.getAttribute('how/stopepochs')
    else:
        start_epoch = time.mktime(datetime.datetime.strptime(obj.startdate+obj.starttime, "%Y%m%d%H%M%S").timetuple())
        end_epoch = time.mktime(datetime.datetime.strptime(obj.enddate+obj.endtime, "%Y%m%d%H%M%S").timetuple())
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
            # This tactic assumes that all quantities are available in all scans
            # of the pvol, and that all quantities are ordered the same way in
            # each scan. This is probably dangerous.
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
                dtype = hfile.getScan(0).getParameter(field_name).getData().dtype
            field_dic['data'] = np.zeros((total_rays, nbins), dtype=dtype)
            start = 0
            # loop over the sweeps, copy data into correct location in data array
            if rio.objectType is _rave.Rave_ObjectType_SCAN:
                sweep_data = obj.getParameter(field_name).getData()
                field_dic['data'][start:start + obj.nrays] = sweep_data[:]
            elif rio.objectType is _rave.Rave_ObjectType_PVOL:
                for i in range(hfile.getNumberOfScans()):
                    scan = hfile.getScan(i)
                    sweep_data = scan.getParameter(field_name).getData()
                    param = scan.getParameter(field_name)
                    scaled_sweep_data = sweep_data[:] * param.gain + param.offset
                    field_dic['data'][start:start + scan.nrays] = scaled_sweep_data[:]
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
