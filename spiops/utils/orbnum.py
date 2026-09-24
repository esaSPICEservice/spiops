from spiops.utils.files import list_files_from_ftp, download_file
from datetime import datetime
import os


def get_orbnum(mission, pattern):
    ftp_orbnum_folder = f'/data/SPICE/{mission}/misc/orbnum/'
    files = list_files_from_ftp(ftp_orbnum_folder, pattern)
    files.sort()
    print(files)
    if len(files) > 0:
        file = files[-1]
        if os.path.isfile(file):
            print('File already downloaded: ' + file)
        else:
            download_file(ftp_orbnum_folder, file)
        return file
    return None

def get_orbnums(mission, pattern):
    ftp_orbnum_folder = f'/data/SPICE/{mission}/misc/orbnum/'
    files = list_files_from_ftp(ftp_orbnum_folder, pattern)
    files.sort()
    if len(files) > 0:
        for file in files:
            if os.path.isfile(file):
                print('File already downloaded: ' + file)
                continue
            download_file(ftp_orbnum_folder, file)
        return files
    return None


def parse_bc_orbnum(orbnum_path):

    headers = [
    "No.", "Event_UTC_APO", "Event_SCLK_APO", "OP-Event_UTC_PERI", 
    "SolLon", "SolLat", "SC_Lon", "SC_Lat", "Alt", "Inc", "Ecc", 
    "LonNode", "Arg_Per", "Sol_Dist", "Semi_Axis"
    ]
    data_rows = []
    with open(orbnum_path, mode="r", encoding="utf-8") as file:
        # Read all lines from the file
        lines = file.readlines()
        
        # Skip the first two rows (original header text and the ===== separator line)
        for line in lines[2:]:
            # Split by any whitespace and strip trailing newlines
            row_data = line.split()
            
            if row_data:
                # Convert numeric strings to float/int where applicable
                processed_row = []
                for item in row_data:
                    try:
                        # Convert to integer if it has no decimal point
                        if '.' not in item:
                            processed_row.append(int(item))
                        else:
                            processed_row.append(float(item))
                    except ValueError:
                        # Keep as string if it is a timestamp or SCLK string
                        processed_row.append(item)
                
                data_rows.append(processed_row)
    return headers, data_rows

class BCOrbnumHandler:

    def __init__(self, orbnum_path):
        self.body_radius = 2439.7
        self.headers, self.raw_data = parse_bc_orbnum(orbnum_path)
        self._extend_drifts("OP-Event_UTC_PERI", "PERI_Drift")
        self._extend_drifts("Event_UTC_APO", "APO_Drift")
        self._extend_periapsis_altitude()
        self.filtered_data = list(self.raw_data)
        
    def _get_col_index(self, header_name):
        """Helper to fetch column index from header string."""
        return self.headers.index(header_name)

    def _get_column_values(self, header_name):
        """Helper to extract a specific column from the active filtered dataset."""
        idx = self._get_col_index(header_name)
        return [row[idx] for row in self.filtered_data]

    def _to_dt(self, t_str):
        """Helper to parse the file timestamp into a comparable datetime object."""
        return datetime.fromisoformat(t_str.replace('Z', ''))

    def _extend_drifts(self, event, col_name):
        event_index = self._get_col_index(event)
        event_data = []
        for row in self.raw_data:
            event_data.append(row[event_index])
        timestamps = [datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ") for t in event_data]
        drifts = []
        diff = 0
        for i in range(len(timestamps) - 1):
            diff = (timestamps[i+1] - timestamps[i]).total_seconds()
            drifts.append(diff)
        # Add the last difference at the end of the drift array
        drifts.append(diff)
        for row, drift in zip(self.raw_data, drifts):
            row.append(drift)
        self.headers.extend([col_name])

    def _extend_periapsis_altitude(self):
        peri_altitudes = []
        for row in self.raw_data:
            peri_altitude = self._calculate_peripasis_altitude(
                row[self._get_col_index("Ecc")],
                row[self._get_col_index("Semi_Axis")])
            peri_altitudes.append(peri_altitude)

        for row, drift in zip(self.raw_data, peri_altitudes):
            row.append(drift)
        self.headers.extend(["PERI_Alt"])

    def _calculate_peripasis_altitude(self, ecc, a):
        return (a * (1 - ecc)) - self.body_radius

    # ==========================================
    # State Management Methods
    # ==========================================
    def reset(self):
        """Resets the filter state back to the original dataset."""
        self.filtered_data = list(self.raw_data)
        return self

    # ==========================================
    # Filter Methods (Chainable)
    # ==========================================
    def min_orbit(self, val):
        """Filters rows where orbit number ('No.') >= val."""
        idx = self._get_col_index("No.")
        self.filtered_data = [row for row in self.filtered_data if row[idx] >= val]
        return self

    def max_orbit(self, val):
        """Filters rows where orbit number ('No.') <= val."""
        idx = self._get_col_index("No.")
        self.filtered_data = [row for row in self.filtered_data if row[idx] <= val]
        return self

    def min_date(self, start_iso):
        """Filters rows where Event_UTC_APO >= start_iso string (e.g., '2027-03-14T05:00:00Z')."""
        idx = self._get_col_index("Event_UTC_APO")
        start_dt = self._to_dt(start_iso)
        self.filtered_data = [row for row in self.filtered_data if self._to_dt(row[idx]) >= start_dt]
        return self

    def max_date(self, end_iso):
        """Filters rows where Event_UTC_APO <= end_iso string (e.g., '2027-03-14T18:00:00Z')."""
        idx = self._get_col_index("Event_UTC_APO")
        end_dt = self._to_dt(end_iso)
        self.filtered_data = [row for row in self.filtered_data if self._to_dt(row[idx]) <= end_dt]
        return self

    # ==========================================
    # Accessor / Terminal Methods
    # ==========================================
    def get_event_utc_apo(self):
        return self._get_column_values("Event_UTC_APO")

    def get_op_event_utc_peri(self):
        return self._get_column_values("OP-Event_UTC_PERI")

    def get_apo_alt(self):
        return self._get_column_values("Alt")

    def get_inclination(self):
        return self._get_column_values("Inc")

    def get_ecc(self):
        return self._get_column_values("Ecc")

    def get_semi_axis(self):
        return self._get_column_values("Semi_Axis")

    def get_lon_node(self):
        return self._get_column_values("LonNode")

    def get_arg_per(self):
        return self._get_column_values("Arg_Per")
    
    def get_peri_drift(self):
        return self._get_column_values("PERI_Drift")
    
    def get_apo_drift(self):
        return self._get_column_values("APO_Drift")
    
    def get_peri_altitude(self):
        return self._get_column_values("PERI_Alt")
    
    def get_orbit_number(self):
        return self._get_column_values("No.")

    def get_data(self):
        """Returns the full row structures matching current filters."""
        return self.filtered_data

    def get_orbital_parameters(self, drifts=False):

        if drifts:
            return {
                'Periapsis drift': {
                    'data': self.get_peri_drift(),
                    'units': 'seconds'
                },
                'Apoapsis drift': {
                    'data': self.get_apo_drift(),
                    'units': 'seconds'
                }
            }



        return {
        'Apoapsis altitude': {
            'data': self.get_apo_alt(),
            'units': 'km'
        },
        'Periapsis altitude': {
            'data': self.get_peri_altitude(),
            'units': 'km'
        },
        'Inclination': {
            'data': self.get_inclination(),
            'units': 'degrees'
        },
        'Eccentricity': {
            'data': self.get_ecc(),
            'units': ''
        },
        'Semi-major axis': {
            'data': self.get_semi_axis(),
            'units': 'km'
        },
        'Longitude of ascending node': {
            'data': self.get_lon_node(),
            'units': 'degrees'
        },
        'Argument of periapsis': {
            'data': self.get_arg_per(),
            'units': 'degrees'
        }
    }