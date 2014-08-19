
import _raveio
import matplotlib.pyplot as plt

def plot(data, ax):
    ax.set_title("B-scan")
    im = ax.imshow(data, cmap='jet', clim=(0,255))
    plt.colorbar(im)

# read and extract data from file
rio = _raveio.open("data/Example_pvol.h5")
pvol = rio.object
scan = pvol.getScan(0)
dbzh = scan.getParameter("DBZH")
dbzh_data = dbzh.getData()

# plot
fig = plt.figure(figsize=(16,12))
ax = fig.add_subplot(211)
plot(dbzh_data, ax=ax)

#############################################################
# read using RAVE, convert to Py-ART, convert to RAVE, plot #
#############################################################

# read using RAVE
import _raveio
rio2 = _raveio.open("data/Example_pvol.h5")

# convert to Py-ART Radar object
import rave_pyart
radar = rave_pyart.toPyART(rio2)

# convert to RAVE object
import pyart_rave
rio3 = pyart_rave.fromPyART(radar)

# plot
ax = fig.add_subplot(212)
pvol3 = rio3.object
#pvol3 = rio.object
scan3 = pvol3.getScan(0)
dbzh3 = scan3.getParameter("DBZH")
dbzh_data3 = dbzh3.getData()
plot(dbzh_data3, ax=ax)


plt.show()
