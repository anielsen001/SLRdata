from __future__ import print_function

import sys
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

import numpy as np
import pandas as pd
from datetime import datetime

# CRD file parser

# Header lines H1 to H9 define one "unit"
# Header lines H4 to H8 define one session

# station metadata from H2
# target metadata from H3
# Config metadata C0, C1, C2, C3, C4
# All of these metadata can be either in session or unit.
# Have to ask session for metadata, and if not available, delegate request to unit.
# TODO: Load metadata in order, warn if session defines metadata already in unit.

# Inside session, time series data comes from record types 10, 11, and 12 (ranges), 
# 20 and 21 (meteo), 30 (angles) and 40 (calibration).
# Also type 50 (pass statistics) can appear.
# The time series stuff should go in a Pandas dataframe?

# References:
# https://edc.dgfi.tum.de/en/oc/crd/
# https://ilrs.cddis.eosdis.nasa.gov/docs/2009/crd_v1.01.pdf

def _int( x ):
    """
    string in this code are cast using "int", but in come cases a
    valid string may be left empty and no number is required. In this
    case, return a nan, otherwise raise the error
    """
    try:
        return int( x )
    except ValueError:
        if len( x.strip() ) > 0 : raise
        else: return np.nan

config_types = { '0': 'system',
                 '1': 'laser',
                 '2': 'detector',
                 '3': 'timing',
                 '4': 'transponder' }

def parse_unit(line):
    """Parse a H1 line into a unit dict."""
    if not line.startswith("H1"):
        raise ValueError("Not H1 line for unit parser!")
    timestamp = {
        "year" : int(line[10:14]),
        "month" : int(line[15:17]),
        "day" : int(line[18:20]),
        "hour" : int(line[21:23])
    }
    return {
        "format" : line[3:6],
        "version" : line[7:9].strip(),
        "time" : datetime(**timestamp),
        "sessions" : []
    }

def parse_session(line):
    """Parse H4 line into session."""
    if not line.startswith("H4"):
        raise ValueError("Not H4 line for session parser!")
    res = {
        "data" : None,
        "start" : datetime(
            year = int(line[6:10]),
            month = int(line[11:13]),
            day = int(line[14:16]),
            hour = int(line[17:19]),
            minute = int(line[20:22]),
            second = int(line[23:25])
        ),
        "troposphere_corrected" :  int(line[49]),
        "CoM_corrected" : int(line[51]),
        "receive_amplitude_corrected" : int(line[53]),
        "station_delay_corrected" : int(line[55]),
        "spacecraft_delay_corrected" : int(line[57]),
        "range_type" : int(line[59]),
        "data_quality" : int(line[61])
    }
    try:
        # end time is not required by CRD format. If not set
        # all values should be set to -1 
        res["end"] = datetime(
            year = int(line[26:30]),
            month = int(line[31:33]),
            day = int(line[34:36]),
            hour = int(line[37:39]),
            minute = int(line[40:42]),
            second = int(line[43:45])
        )
    except ValueError:
        # use NaN to represent not available
        res['end'] = np.nan
    return res


def parse_station(line):
    if not line.startswith("H2"):
        raise ValueError("Not H2 line for Station constructor!")
    return {
        "name" : line[3:13].strip(),
        "ID" : int(line[14:18]),
        "system" : int(line[19:21]),
        "occupancy" : int(line[22:24]),
        "timescale" : int(line[25:27])
    }


def parse_target(line):
    """Parse H3 line for target parameters."""
    if not line.startswith("H3"):
        raise ValueError("Not H3 line for Target constructor!")
    return {
        "name" : line[3:13].strip(),
        "ID" : int(line[14:22]),
        "SIC" : int(line[23:27]),
        "NORAD" : _int(line[28:36]),
        "timescale" : int(line[37:38]),
        "type" : int(line[39])
    }

# the configuration lines are free format. White space parsing
# is required
def parse_system_configuration( line ):
    """
    parse the C0 system configuration line
    these lines may be repeated
    """
    
    if not line.startswith("C0"):
        raise ValueError('Not C0 line for configuration')
    
    # parse the line by splitting
    parts = line.split()

    record_type = parts[0]
    detail_type = int( parts[1] )
    transmit_wavelength_nm = float( parts[2] )

    system_config_id = parts[3]

    # these may not all be present
    try:
        component_a_config_id = parts[4]
    except IndexError:
        component_a_config_id = np.nan
    try:
        component_b_config_id = parts[5]
    except IndexError:
        component_b_config_id = np.nan
    try:
        component_c_config_id = parts[6]
    except IndexError:
        component_c_config_id = np.nan
    try:
        component_d_config_id = parts[7]
    except IndexError:
        component_d_config_id = np.nan
                
    res = { 'record_type' : record_type,
             'detail_type' : detail_type,
             'transmit_wavelength_nm' : transmit_wavelength_nm,
             'system_config_id' : system_config_id,
             'component_a_config_id' : component_a_config_id,
             'component_b_config_id' : component_b_config_id,
             'component_c_config_id' : component_c_config_id,
             'component_d_config_id' : component_d_config_id
             }
    
    return res

def parse_laser_configuration( line ):
    """
    parse the laser configuration line

    there may be more than one line if e.g. two-color ranging
    is being used
    """

    if not line.startswith("C1"):
        raise ValueError('Not C1 line for configuration')
    
    parts = line.split()

    record_type = parts[0]
    detail_type = int( parts[1] )
    laser_config_id = parts[2]
    laser_type = parts[3]
    primary_wavelength_nm = float( parts[4] )
    nominal_fire_rate_hz = float( parts[5] )
    pulse_energy_mj = float( parts[6] )       # should be updated when changes by 10%
    pulse_width_ps = float( parts[7] )        # should be updated when changes by 10%
    beam_divergence_arcsec = float( parts[8] )
    num_pulses_in_semitrain = int( parts[9] )

    res = { 'record_type' : record_type,
            'detail_type' : detail_type,
            'laser_config_id' : laser_config_id,
            'laser_type' : laser_type,
            'primary_wavelength_nm' : primary_wavelength_nm,
            'nominal_fire_rate_hz' : nominal_fire_rate_hz,
            'pulse_energy_mj' : pulse_energy_mj,
            'pulse_width_ps' : pulse_width_ps,
            'beam_divergence_arcsec' : beam_divergence_arcsec,
            'num_pulses_in_semitrain' : num_pulses_in_semitrain
            }
    return res

def parse_detector_configuration( line ):
    if not line.startswith("C2"):
        raise ValueError('Not C2 line for configuration')
    
    parts = line.split()

    record_type = parts[0]
    detail_type = int( parts[1] )
    detector_config_id = parts[2]
    detector_type = parts[3]
    applicable_wavelength_nm = float( parts[4] )
    quantum_efficiency = float( parts[5] )
    applied_voltage = np.nan if parts[6].lower() == 'unknown' else float( parts[6] )
    dark_count_khz = np.nan if parts[7].lower() == 'unknown' else float( parts[7] )
    output_pulse_type = parts[8]
    output_pulse_width = np.nan if parts[9].lower() == 'unknown' else float( parts[9] )
    spectral_filter_nm = float( parts[10] )
    trans_percent_spectral_filter = float( parts[11] )
    spatial_filter_arcsec = np.nan if parts[12].lower() == 'unknown' else float( parts[12] )
    ext_signal_proc = parts[13]

    res = { 'record_type' : record_type,
            'detail_type' : detail_type,
            'detector_config_id' : detector_config_id,
            'detector_type' : detector_type,
            'applicable_wavelength_nm' : applicable_wavelength_nm,
            'quantum_efficiency' : quantum_efficiency,
            'applied_voltage' : applied_voltage,
            'dark_count_khz' : dark_count_khz,
            'output_pulse_type' : output_pulse_type,
            'output_pulse_width' : output_pulse_width,
            'spectral_filter_nm' : spectral_filter_nm,
            'trans_percent_spectral_filter' : trans_percent_spectral_filter,
            'spatial_filter_arcsec' : spatial_filter_arcsec,
            'ext_signal_proc' : ext_signal_proc
            }
    
    return res

def parse_timing_configuration( line ):
    if not line.startswith("C3"):
        raise ValueError('Not C3 line for configuration')
    
    parts = line.split()

    record_type = parts[0]
    detail_type = int( parts[1] )
    timing_config_id = parts[2]
    time_source = parts[3]
    frequency_source = parts[4]
    timer = parts[5]
    timer_serial_number = parts[6]
    epoch_delay_correction_us = float( parts[7] )

    res = {
        'record_type' : record_type,
        'detail_type' : detail_type,
        'timing_config_id' : timing_config_id,
        'time_source' : time_source,
        'frequency_source' : frequency_source,
        'timer' : timer,
        'timer_serial_number' : timer_serial_number,
        'epoch_delay_correction_us' : epoch_delay_correction_us
        }
    
    return res

def parse_transponder_configuration( line ):
    if not line.startswith("C4"):
        raise ValueError('Not C4 line for configuration')
    
    parts = line.split()
 
    record_type = parts[0]
    detail_type = int( parts[1] )
    transponder_config_id = parts[2]
    estimated_station_utc_offset_ns = float( parts[3] )
    estimated_station_osc_drift = float( parts[4] )
    estimated_transponder_utc_offset_ns = float( parts[5] )
    estimated_transponder_osc_drift = float( parts[6] )
    transponder_clock_ref_time = float( parts[7] )
    station_clock_correction_applied = int( parts[8] )
    spacecraft_clock_correction_applied = int( parts[9] )
    spacecraft_time_simplified = int( parts[10] )

    res = {
        'record_type' : record_type,
        'detail_type' : detail_type,
        'transponder_config_id' : transponder_config_id,
        'estimated_station_utc_offset_ns' : estimated_station_utc_offset_ns,
        'estimated_station_osc_drift' : estimated_station_osc_drift,
        'estimated_transponder_utc_offset_ns' : estimated_transponder_utc_offset_ns,
        'estimated_transponder_osc_drift' : estimated_transponder_osc_drift,
        'transponder_clock_ref_time' : transponder_clock_ref_time,
        'station_clock_correction_applied' : station_clock_correction_applied,
        'spacecraft_clock_correction_applied' : spacecraft_clock_correction_applied,
        'spacecraft_time_simplified' : spacecraft_time_simplified
        }
    return res

config_parse_funcs = { 'laser' : parse_laser_configuration,
                       'detector' : parse_detector_configuration,
                       'timing' : parse_timing_configuration,
                       'transponder' : parse_transponder_configuration }

class ConfigRecord( object ):
    """
    handles the config records as a dictionary of lists of dictionaries
    """

    # configuration id
    system_config_id = None

    # laser transmit wavelength, nm
    laser_wavelength_nm = None

    # these are the ids of the other components
    # created dynamically at startup based on _comp_types
    #laser_config_id = None
    #detector_config_id = None
    #timing_config_id = None
    #transponder_config_id = None

    # The documentation is not clear that any of these map to any particular
    # component, so set these on init and as subsequent config lines are match
    # set those ids
    _comp_id = {'a': None,
                'b': None,
                'c': None,
                'd': None }

    # these are the component types
    _comp_types = [ 'laser',
                    'detector',
                    'timing',
                    'transponder' ]

    # keep the line that was parsed for referenced
    _line = None

    def __init__( self, *args, **kwargs ):
        # create component type id class variables
        for comp in self._comp_types:
            setattr( self, comp + '_config_id', None )
            setattr( self, comp + '_config', None )
        
        # initialize with c0 line from crd file
        self._line = args[0]

        rec = parse_system_configuration( self._line )
        self.system_config_id = rec[ 'system_config_id' ]

        self.laser_wavelength_nm = rec[ 'transmit_wavelength_nm' ]

        self._comp_id['a'] = rec[ 'component_a_config_id' ]
        self._comp_id['b'] = rec[ 'component_b_config_id' ]
        self._comp_id['c'] = rec[ 'component_c_config_id' ]
        self._comp_id['d'] = rec[ 'component_d_config_id' ]
        
        # no return value

    def add_component( self, component, component_config ):
        """
        add a component configuration to the system configuration
        """
        # construct th compenent id method variable name
        comp_id_var = component + '_config_id'
        config_var = component + '_config'
        
        # check if this config already has this component config
        if getattr( self, comp_id_var ) is not None:
            raise ValueError( str( self.system_config_id) + \
                              'already has a ' + component + ' config (' +\
                              getattr( selfk, comp_id_var) + ')' )

        # check if the component_config_id matches any of the component ids
        match = False
        for k,v in self._comp_id.items():
            if v.lower() in component_config[ comp_id_var ].lower():
                # good match, set values
                match = True
                setattr( self, comp_id_var, component_config[ comp_id_var ] )
                setattr( self, config_var, component_config )
                break

        if not match:
            raise ValueError( component_config[ comp_id_var ] +
                              ' did not match any components in ' +
                              self.system_config_id )

class DataError(Exception):
    pass
        
class Data(object):
    """
    Data is a generic class to hold sampled data. 
    This will include the engineering, normal point and meterological 
    data

    The object should be used in a read mode, parse mode and then
    access mode.
    """
    # line start code
    line_start_code = ''

    # this is a list of the field names for this data type
    field_names = []
    
    # this will be a large string of text parsed from each
    # record in the file, until they are all accumulated and
    # ready to parse into a data frame
    _lines = ''

    # load each line, after loading each line it fill be put
    # into a pandas data frame
    _df = None

    def __init__(self):

        # create the first line in _lines to be the header
        header_line = ' '.join( self.field_names ) + '\n'
        self._lines += header_line + '\n'

    def add_line(self,line):
        """ get the line of data add concatenate to _lines """
        # on read, the lines are split with '\n' so add it back
        # here
        if not line.startswith(self.line_start_code):
            raise DataError('wrong line type passed' )
        self._lines += line + '\n'

    def parse(self):
        """ parse the collected data into a pandas data frame """
        self._df = pd.read_csv( StringIO(self._lines), delim_whitespace=True )

    def __getitem__( self, item ):
        """ extract collected data """
        return self._df[item]

        
class NormalPointData(Data):
    line_start_code = '11'    
    field_names = ['record_type',
                   'seconds_of_day',
                   'time_of_flight',
                   'system_config_id',
                   'epoch_event',
                   'nomal_point_window_length',
                   'number_raw_ranges',
                   'bin_rms',
                   'bin_skew',
                   'bin_kurtosis',
                   'bin_peak',
                   'return_rate',
                   'detector_channel']

        
class MeteorologicalData(Data):
    line_start_code = '20'
    field_names = ['record_type',
                   'seconds_of_day',
                   'surface_pressure_mbar',
                   'surface_temperature_kelvin',
                   'relative_humidity_percent',
                   'origin_of_values' ]

class MeteorologicalSupplementData(Data):
    line_start_code = '21'
    field_names = ['record_type',
                   'seconds_of_day',
                   'wind_speed_mps',
                   'wind_direction_az_deg',
                   'precipitation_type',
                   'visibility_km',
                   'sky_clarity',
                   'atmospheric_seeing',
                   'cloud_cover_percent' ]
                   
    
        
def parse_CRD(data):
    active_unit = None
    active_session = None
    active_data = []
    active_npt_data = None
    active_met_data = None
    active_met_supp_data = None
    units = []

    # track the configuration status
    # put a flag to track if the config belongs to the unit
    # or the session
    unit_config = False
    session_config = False
    parsing_config = False  # assume all the config lines are sequential (??)
    active_conf = None      # this is the active configuraiton object
    
    for line in data.split("\n"):
        line = line.upper()

        # if we have been parsing a config section, need to determine
        # if we are done
        if parsing_config and not line.startswith('C'):
            # assume done parsing the config
            parsing_config = False
            
            # associate with session or unit
            if not active_session:
                active_unit['config'] = active_conf
            else:
                active_session['config'] = active_conf

            active_conf = None

        # the H lines are header line
        if line.startswith("H1"):
            # Start of unit.
            active_unit = parse_unit(line)
            units.append(active_unit)
            parsing_config = False
            
        elif line.startswith("H9"):
            # End of unit.
            active_unit = None
            parsing_config = False
            
        elif line.startswith("H4"):
            # Start of session, add new session to active unit.
            active_session = parse_session(line)
            active_unit["sessions"].append(active_session)
            parsing_config = False
            
        elif line.startswith("H8"):
            # End of session, convert active_data into array and save to
            # active session.
            active_session["data"] = np.array(active_data)
            active_data = []
            
            active_session['npt'] = active_npt_data
            if active_session['npt']:
                active_session['npt'].parse()
            active_npt_data = None
    
            active_session['met'] = active_met_data
            if active_session['met']:
                active_session['met'].parse()
            active_met_data = None
            
            active_session['met_supp'] = active_met_supp_data
            if active_session['met_supp']:
                active_session['met_supp'].parse()
            active_met_supp_data = None

            active_session = None
            parsing_config = False
            
        elif line.startswith("H2"):
            # Station definition, add station to active session, 
            # or to active unit if no active session.
            if active_session is None:
                active_unit["station"] = parse_station(line)
            else:
                active_session["station"] = parse_station(line)
            parsing_config = False
                
        elif line.startswith("H3"):
            # Target definition, add new target to active session, 
            # or to active unit if no active session.
            if active_session is None:
                active_unit["target"] = parse_target(line)
            else:
                active_session["target"] = parse_target(line)
            parsing_config = False

        # configuration records can belong to either the session
        # or the unit
        elif line.startswith("C0"):
            # configuration line 0
            system_config = ConfigRecord( line )

            parsing_config = True
            
            # there may be more than one configuration used at a time
            # if multiple lasers are used simultaneously
            if not active_conf:
                active_conf = [ system_config ] # put in list 
            else:
                active_conf.append( system_config )

        elif line.startswith('C'):
            # already caught C0, so don't worry about it.
            if not parsing_config:
                raise ValueError('Not actively parsing config')

            linetype = config_types[ line[1] ]
            parse_func = config_parse_funcs[ linetype ]
            comp_config = parse_func( line )
            # put into the system config
            # the system config may be a list of length > 1
            # this could match potentially all of the configurations
            for ac in active_conf:
                ac.add_component( linetype, comp_config )            
            
        elif line.startswith("10"):
            # full-rate/engineering/quicklook data
            # Data point, add to active_data.
            sline = line.split()
            t = float(sline[1])
            r = float(sline[2])
            active_data.append([t, r])
            parsing_config = False

        elif line.startswith('11'):
            # normal point data
            if not active_npt_data:
                active_npt_data = NormalPointData()
            active_npt_data.add_line(line)
            parsing_config = False
            
        elif line.startswith('20'):
            if not active_met_data:
               active_met_data = MeteorologicalData()
            active_met_data.add_line(line)            
            parsing_config = False

        elif line.startswith('21'):
            if not active_met_supp_data:
                active_met_supp_data = MeteorologicalSupplementData()
            active_met_supp_data.add_line(line)            
            parsing_config = False
            
        else:
            # there are other types of lines that are not handled,
            # handle those if we encounter them.
            parsing_config = False
            
    return units
        




if __name__=="__main__":
    from sys import argv
    import json
    with open(argv[1]) as f:
        Units = parse_CRD(f.read())
        print(json.dumps(Units[0], default=str, indent=4))

