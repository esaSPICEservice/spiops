import spiceypy as spice
from spiceypy.utils.support_types import SpiceyError
from datetime import timezone


def get_timeline_events(sc):

    try:
        """
        Load the JUICE timeline kernel and return a list or dict of event info.
        """
        # 1. Retrieve the string arrays from the kernel
        #    spice.gcpool returns a Python list of strings
        event_types = spice.gcpool(sc + "_TIMELINE_EVENT_TYPE", 0, 9999)
        event_names = spice.gcpool(sc + "_TIMELINE_EVENT_NAME", 0, 9999)
        event_times = spice.gcpool(sc + "_TIMELINE_EVENT_TIME", 0, 9999)

        # 2. Create a data structure. For example, a list of dictionaries
        events = []
        for etype, ename, etime in zip(event_types, event_names, event_times):
            # Convert time string --> ephemeris time (ET)
            event_et = spice.str2et(etime)
            # Convert ephemeris time --> Python datetime in UTC
            event_utc = spice.et2datetime(event_et)

            events.append({
                "type": etype,
                "name": ename,
                "datetime_utc": event_utc,
                "ephemeris_time": event_et
            })

        # 3. Sort events by naive datetime
        events = sorted(events, key=lambda e: e["datetime_utc"])

        return events

    except SpiceyError:
        print('Error retrieving ' + sc + " timeline events.")
        return None


def get_closest_event(t, events, datetime_key="datetime_utc"):
    min_diff = None
    closest_event = None
    for ev in events:
        # Also ensure the event times are naive
        ev_time = ev[datetime_key]
        event_utc = ev_time.astimezone(timezone.utc)

        if ev_time.tzinfo is not None:
            ev_time = ev_time.replace(tzinfo=None)

        diff_seconds = abs((ev_time - t).total_seconds())
        if min_diff is None or diff_seconds < min_diff:
            min_diff = diff_seconds
            closest_event = ev

    return closest_event
