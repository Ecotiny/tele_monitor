from nclient import Client, dec_to_dms
import numpy as np
import PyIndi
import time

class IndiClient(PyIndi.BaseClient):
    def __init__(self):
        super(IndiClient, self).__init__()
    def newDevice(self, d):
        pass
    def newProperty(self, p):
        print(p.getName())
    def removeProperty(self, p):
        pass
    def newBLOB(self, bp):
        pass
    def newSwitch(self, svp):
        pass
    def newNumber(self, nvp):
        pass
    def newText(self, tvp):
        pass
    def newLight(self, lvp):
        pass
    def newMessage(self, d, m):
        pass
    def serverConnected(self):
        pass
    def serverDisconnected(self, code):
        pass
 

if __name__ == "__main__":
    n = 500
    ras = np.linspace(0, 360, n)
    decs = np.linspace(-90, 90, n)
    indiclient = IndiClient()
    indiclient.setServer("localhost",7624)
     
    indiclient.connectServer()
    # connect the scope
    telescope="Celestron GPS"
    device_telescope=None
    telescope_connect=None
     
    # get the telescope device
    device_telescope=indiclient.getDevice(telescope)
    exit()
    while not(device_telescope):
        time.sleep(0.5)
        device_telescope=indiclient.getDevice(telescope)
         
    # wait CONNECTION property be defined for telescope
    telescope_connect=device_telescope.getSwitch("CONNECTION")
    while not(telescope_connect):
        time.sleep(0.5)
        telescope_connect=device_telescope.getSwitch("CONNECTION")
     
    # if the telescope device is not connected, we do connect it
    if not(device_telescope.isConnected()):
        # Property vectors are mapped to iterable Python objects
        # Hence we can access each element of the vector using Python indexing
        # each element of the "CONNECTION" vector is a ISwitch
        telescope_connect[0].s=PyIndi.ISS_ON  # the "CONNECT" switch
        telescope_connect[1].s=PyIndi.ISS_OFF # the "DISCONNECT" switch
        indiclient.sendNewSwitch(telescope_connect) # send this new value to the device

    while not device_telescope.getNumber("EQUATORIAL_EOD_COORD"):
        time.sleep(0.5)

    print("Connected!\nInitialising interface...")
    cli = Client()

    while (1):
        radec = device_telescope.getNumber("EQUATORIAL_EOD_COORD")
        track = device_telescope.getSwitch("TELESCOPE_TRACK_MODE")
        track_type = -1
        for idx, switch in enumerate(track):
            if switch.s == 1: track_type = idx
            
        try:
            if radec:
                # convert hours to degrees
                ra = (radec[0].value * 360)/24
                dec = radec[1].value

                tracking = ['Sidereal', 'Solar', 'Lunar', 'Custom'][track_type]
                rad, ram, rasec = dec_to_dms(ra)
                decd, decm, decsec = dec_to_dms(dec)
                cli.update(ra, dec, {
                    "   RA": f"{rad: 3d} {ram: 2d} {rasec: 4.2f}",
                    "  Dec": f"{decd: 3d} {decm: 2d} {decsec: 4.2f}",
#                    "  Alt": f"{alt:.3f}",
#                    "   Az": f"{az:.3f}",
                    "Scope": telescope,
                    "Track": tracking})
                time.sleep(0.1)
            else:
                raise KeyboardInterrupt
        except Exception as e:
            cli.shutdown()
            print("Exiting...")
            raise e
            break

