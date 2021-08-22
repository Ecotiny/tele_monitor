import curses
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
import astropy.units as u
from astropy.time import Time
from astroquery.simbad import Simbad
Simbad.add_votable_fields("flux(V)")
import numpy as np
import time
import json
import os
import sqlite3

def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate

def ra_dec_to_map(scrn, ra, dec):
    rows, cols = scrn.getmaxyx()
    rows += -2 # to allow for the border
    cols += -2

    center = (0, 0) # deg east, deg south
    ra = (ra + center[0])%360
    dec = (((dec + center[1])+90) % 180) - 90

    # ra will be from 0 to 360
    # dec will be from -90 -> +90
    # ra will be on the x axis, dec on the y
    rapx = ((ra/360)*cols) + 1
    decpx = (((dec + 90) / 180)*rows) + 1

    return int(rapx), int(decpx)

def dec_to_dms(dec):
    degrees = dec // 1
    dec -= degrees
    dec *= 60
    minutes = dec // 1
    dec -= minutes
    dec *= 60
    seconds = dec
    return (int(degrees), int(minutes), seconds)

def dms_to_dec(dms):
    if type(dms) == str:
        dms = [*map(float, dms.strip().split(" "))]
    d, m, s = dms
    dec = d
    dec += m / 60
    dec += s / 3600
    return dec 

def hms_to_dec(hms):
    dec = dms_to_dec(hms) * 360/24
    return dec 

def dec_to_hms(dec):
    return dec_to_dms(dec*24/360)

# make this static, only update when time is significantly delta if it takes too long
@static_vars(horizon=None, updated=Time.now())
def draw_horizon(scrn):
    if not draw_horizon.horizon or Time.now() - draw_horizon.updated  > 1 * u.minute:
        location = EarthLocation(lon=-170.2 * u.deg, lat=-45.5 * u.deg, height=100 * u.m)
        obstime = Time.now() 
        draw_horizon.updated = obstime
        n = 100
        az_list = np.linspace(0, 360, n)
        altaz_frame = AltAz(location=location, obstime=obstime)
        draw_horizon.horizon = []
        for az in az_list:
            coord = SkyCoord(alt=0 * u.deg, az=az * u.deg, frame=altaz_frame)
    
            radec = coord.transform_to('icrs')
            ra = radec.ra.deg
            dec = radec.dec.deg
            rapx, decpx = ra_dec_to_map(scrn, ra, dec)
            draw_horizon.horizon.append((rapx, decpx))
    for ra, dec in draw_horizon.horizon:
        scrn.addch(dec, ra, "-", curses.color_pair(3))    

def draw_labels(scrn):
    rows, cols = scrn.getmaxyx()
    scrn.addstr(rows-1, cols//2-1, "  N  ")
    scrn.addstr(0, cols//2-1,      "  S  ")
    scrn.addstr(rows//2-1, cols-1, " ") 
    scrn.addstr(rows//2, cols-1,   "W") 
    scrn.addstr(rows//2+1, cols-1, " ") 
    scrn.addstr(rows//2-1, 0,      " ")
    scrn.addstr(rows//2, 0,        "E")
    scrn.addstr(rows//2+1, 0,      " ")

@static_vars(objs=None, bstars=None)
def draw_objects(scrn, curs, maglimit=5):
    if not draw_objects.bstars:
        outstars = []
        for star in curs.execute(f"SELECT * FROM catalogue WHERE mag < {maglimit};"):
            tyc2, ra, dec, mag = star
            rapx, decpx = ra_dec_to_map(scrn, ra, dec)
            outstars.append({"rapx": rapx, "decpx": decpx, "mag": mag})
        draw_objects.bstars = {}
        draw_objects.bstars['stars'] = outstars

    for star in draw_objects.bstars['stars']:
        mag_chart = ["0", "O", "o", "."]
        char = mag_chart[0]
        if star['mag'] > 0:
            char = mag_chart[1]
        if star['mag'] > 1:
            char = mag_chart[2]
        if star['mag'] > 4:
            char = mag_chart[3]
        scrn.addch(star['decpx'], star['rapx'], char, curses.color_pair(1))

    if not draw_objects.objs:
        with open("objs.json", "r") as f:
            f.seek(0, 2)
            if f.tell() > 5: 
                f.seek(0,0)
                draw_objects.objs = json.load(f)
            else:
                draw_objects.objs = "NULL" 
    else:
        if draw_objects.objs != "NULL":
            for obj in draw_objects.objs:
                coord = SkyCoord(obj['coords'], frame='icrs', unit=(u.hour, u.deg))
                rapx, decpx = ra_dec_to_map(scrn, coord.ra.deg, coord.dec.deg)
                scrn.addch(decpx, rapx, obj['char'], curses.color_pair(1))
        
@static_vars(meridian=None, updated=Time.now())
def draw_meridian(scrn):
    if not draw_meridian.meridian or Time.now() - draw_meridian.updated  > 1 * u.minute:
        location = EarthLocation(lon=-170.491486 * u.deg, lat=-45.872402 * u.deg, height=100 * u.m)
        obstime = Time.now()
        draw_meridian.updated = obstime
        n = 30
        alt_list = np.linspace(0, 90, n)
        altaz_frame = AltAz(location=location, obstime=obstime)
        draw_meridian.meridian = []
        for az in [0, 180]: # due north, due south
            for alt in alt_list:
                coord = SkyCoord(alt=alt * u.deg, az=az * u.deg, frame=altaz_frame)
        
                radec = coord.transform_to('icrs')
                ra = radec.ra.deg
                dec = radec.dec.deg
                rapx, decpx = ra_dec_to_map(scrn, ra, dec)
                draw_meridian.meridian.append((rapx, decpx))
    for ra, dec in draw_meridian.meridian:
        scrn.addch(dec, ra, "|", curses.color_pair(2))    

@static_vars(history=[])
def refresh_starmap(scrn, ra, dec, curs):
    # draw position of telescope
    newrapx, newdecpx = ra_dec_to_map(scrn, ra, dec)
    scrn.clear()
    for rapx, decpx in refresh_starmap.history: 
            scrn.addch(decpx, rapx, "+", curses.color_pair(1)) # color pair 1 is white 

    refresh_starmap.history.append((newrapx, newdecpx))
    # cap length at 10
    if len(refresh_starmap.history) > 10:
        refresh_starmap.history = refresh_starmap.history[1:]

    draw_objects(scrn, curs)
    scrn.border()
    draw_labels(scrn)
    draw_meridian(scrn)
    draw_horizon(scrn)
    scrn.addch(newdecpx, newrapx, "O", curses.color_pair(2)) # color pair 2 is red

def refresh_infopanel(scrn, info):
    # info is a dict containing key value pairs
    scrn.clear()
    scrn.border()
    scrn.addstr(0, 2, " Information ", curses.color_pair(1))
    for pos, field in enumerate(info):
        scrn.addstr(2+pos, 2, field, curses.color_pair(4)) # highlight
        scrn.addstr(2+pos, 3+len(field), info[field], curses.color_pair(1))

@static_vars(stars=[], prev_radec=())
def zoom_map(scrn, ra, dec, width, cursor):
    # width is given in degrees
    rows, cols = scrn.getmaxyx()
    height = (width/cols) * rows
    origin_ra = ra  - width/2
    origin_de = dec# - height/2
    query = f"SELECT * FROM catalogue WHERE mRAdeg > {origin_ra} AND mRAdeg < {origin_ra + width} AND mDEdeg > {origin_de} AND mDEdeg < {origin_de + height} AND mag < 10"
    scrn.clear()
    if zoom_map.prev_radec != (origin_ra, origin_de):
        zoom_map.prev_radec = (origin_ra, origin_de)
        stars = []
        for star in cursor.execute(query):
            tyc2, sra, sdec, mag = star
            d_ra = sra - origin_ra
            d_de = sdec - origin_de
            # scale up the delta values to +- 90 and 0-360
            d_ra = d_ra * 360 / width
            d_de = d_de * 180 / height

            # map to screen
            d_ra, d_de = ra_dec_to_map(scrn, d_ra, d_de)

            char = "0"
            if mag > 4:
                char = "O"
            if mag > 6:
                char = "o"
            if mag > 8:
                char = '.'

            stars.append((d_ra, d_de, char))
        
        zoom_map.stars = stars
   
    for d_ra, d_de, char in zoom_map.stars:
        scrn.addch(d_de, d_ra, char, curses.color_pair(1))

    scrn.border()
    scrn.addstr(0, 2, " A better name", curses.color_pair(1))
    scrn.addch(rows//2, cols//2, "O", curses.color_pair(2))


class Client:
    def __init__(self, sidewidth=40, dbfn='tycho2.db'):
        whole_win = curses.initscr()

        curses.start_color()
        curses.use_default_colors()

        # Star color
        curses.init_pair(1, curses.COLOR_WHITE, curses.COLOR_BLACK)
        
        # Highlight color
        curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
        
        # Line colour
        curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)

        # extra highlighted
        curses.init_pair(4, curses.COLOR_YELLOW, curses.COLOR_BLACK)

        curses.cbreak()
        
        rows, cols = whole_win.getmaxyx()
        self.starmap  = curses.newwin(rows, cols-sidewidth)
        self.infobox  = curses.newwin(rows//2+1, sidewidth, 0, cols-sidewidth)
        self.fieldbox = curses.newwin(rows//2, sidewidth, rows//2+1, cols-sidewidth)
        self.conn = sqlite3.connect(dbfn)
        self.curs = self.conn.cursor()

    def update(self, ra, dec, info):
        refresh_starmap(self.starmap, ra, dec, self.curs)

        rad, ram, rasec = dec_to_dms(ra)
        decd, decm, decsec = dec_to_dms(dec)
        refresh_infopanel(self.infobox, info)
        rah, rahm, rahsec = dec_to_hms(ra)
        zoom_map(self.fieldbox, ra, dec, 5, self.curs)

        self.starmap.refresh()
        self.infobox.refresh()
        self.fieldbox.refresh()

    def shutdown(self):
        self.conn.close()
        curses.endwin()

if __name__ == "__main__":
    n = 500
    ras = np.linspace(0, 360, n)
    decs = np.linspace(-90, 90, n)
    cli = Client()
    try:
        for ra, dec in zip(ras, decs):
            rad, ram, rasec = dec_to_dms(ra)
            decd, decm, decsec = dec_to_dms(dec)
            cli.update(ra, dec, {
                "  RA": f"{rad} {ram} {rasec: 4.2f}",
                " Dec": f"{decd} {decm} {decsec: 4.2f}",
                " Tgt": "Testing object",
                "Trck": "Sidereal"})
    except KeyboardInterrupt:
        cli.shutdown()


