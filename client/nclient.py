import curses
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
import astropy.units as u
from astropy.time import Time
import numpy as np
import time
import json

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
    scrn.addstr(rows//2, cols-1,   "E") 
    scrn.addstr(rows//2+1, cols-1, " ") 
    scrn.addstr(rows//2-1, 0,      " ")
    scrn.addstr(rows//2, 0,        "W")
    scrn.addstr(rows//2+1, 0,      " ")

@static_vars(objs=None, bstars=None)
def draw_objects(scrn):
    if not draw_objects.bstars:
        # read the list and then check
        with open("bstars.json", "r") as f:
            f.seek(0, 2)
            if f.tell() > 5: 
                f.seek(0,0)
                draw_objects.bstars = json.load(f)
            else:
                # get list of 100 brightest stars
                url = "https://raw.githubusercontent.com/aduboisforge/Bright-Star-Catalog-JSON/master/BSC.json"
                # keys are "RA", "DEC", "MAG"
                import requests # we only need it if it's not there already
                print("updating bright star catalogue...")
                r = requests.get(url)
                outstars = []
                for star in r.json():
                    if float(star['MAG']) < 4:
                        coordstr = f"{star['RA']} {star['DEC']}"
                        coord = SkyCoord(coordstr, frame='icrs', unit=(u.hourangle, u.deg))
                        rapx, decpx = ra_dec_to_map(scrn, coord.ra.deg, coord.dec.deg)
                        outstars.append({"coords": coordstr, "rapx": rapx, "decpx": decpx, "mag": float(star['MAG'])})
                draw_objects.bstars = {}
                draw_objects.bstars['stars'] = outstars
                with open("bstars.json", "w") as f:
                    json.dump(draw_objects.bstars, f)

    for star in draw_objects.bstars['stars']:
        if star['mag'] < 4:
            mag_chart = ["0", "O", "o", "."]
            char = mag_chart[0]
            if star['mag'] > 0:
                char = mag_chart[1]
            if star['mag'] > 2:
                char = mag_chart[2]
            if star['mag'] > 3:
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
                coord = SkyCoord(obj['coords'], frame='icrs', unit=(u.deg, u.deg))
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
def refresh_starmap(scrn, ra, dec):
    newrapx, newdecpx = ra_dec_to_map(scrn, ra, dec)
    scrn.clear()
    for rapx, decpx in refresh_starmap.history: 
            scrn.addch(decpx, rapx, "+", curses.color_pair(1)) # color pair 1 is white 

    refresh_starmap.history.append((newrapx, newdecpx))
    # cap length at 10
    if len(refresh_starmap.history) > 10:
        refresh_starmap.history = refresh_starmap.history[1:]

    draw_objects(scrn)
    draw_meridian(scrn)
    scrn.addch(newdecpx, newrapx, "O", curses.color_pair(2)) # color pair 2 is red
    scrn.border()
    draw_labels(scrn)
    draw_horizon(scrn)

    scrn.refresh()

def refresh_infopanel(scrn, info):
    # info is a dict containing:
    # ra, dec, target name, is_exposing, track rates(?)
    scrn.clear()
    scrn.border()
    scrn.addstr(0, 2, " Information ", curses.color_pair(1))
    scrn.addstr(2, 2, f" RA: {info['ra']}", curses.color_pair(1))
    scrn.addstr(3, 2, f"Dec: {info['dec']}", curses.color_pair(1))
    scrn.addstr(5, 2, f"Target: {info['target']}", curses.color_pair(1))
    scrn.refresh()

def dec_to_dms(dec):
    degrees = dec // 1
    dec -= degrees
    dec *= 60
    minutes = dec // 1
    dec -= minutes
    dec *= 60
    seconds = dec
    return (int(degrees), int(minutes), seconds)

if __name__ == "__main__":
    
    whole_win = curses.initscr()
    curses.start_color()
    curses.use_default_colors()

    # Star color
    curses.init_pair(1, curses.COLOR_WHITE, curses.COLOR_BLACK)
    
    # Highlight color
    curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
    
    # Line colour
    curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)

    curses.cbreak()
    
    rows, cols = whole_win.getmaxyx()
    infowidth = 30
    starmap = curses.newwin(rows, (2*cols//3))
    infobox = curses.newwin(rows, cols//3, 0, (2*cols//3))

    n = 500
    ras = np.linspace(0, 360, n)
    decs = np.linspace(-90, 90, n)
    for idx in range(n):
        refresh_starmap(starmap, ras[idx], decs[idx])
        rad, ram, rasec = dec_to_dms(ras[idx])
        decd, decm, decsec = dec_to_dms(decs[idx])
        refresh_infopanel(infobox, {'ra': f"+{rad:02d} {ram:02d}' {rasec:05.2f}\"", 'dec': f"{decd:02d} {decm:02d}' {decsec:05.2f}\"", 'target': idx})
        time.sleep(0.1)
    
