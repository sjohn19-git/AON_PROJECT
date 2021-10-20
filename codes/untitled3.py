import datetime as dt
import matplotlib.dates as mdates
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import dill
from datetime import datetime
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import gridspec
import dill
from matplotlib.dates import MO, TU, WE, TH, FR, SA, SU
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes

fig=plt.figure()

host = fig.add_axes([0.15, 0.1, 0.7, 0.3], axes_class=HostAxes)
par1 = ParasiteAxes(host)

fig.set_figwidth(20)
host.parasites.append(par1)

host.axis["right"].set_visible(False)

par1.axis["right"].set_visible(True)
par1.axis["right"].major_ticklabels.set_visible(True)
par1.axis["right"].label.set_visible(True)

par2.axis["right2"] = par2.new_fixed_axis(loc="right", offset=(60, 0))

p1, = host.plot([0, 1, 2], [0, 1, 2], label="Density")
p2, = par1.plot([0, 1, 2], [0, 3, 2], label="Temperature")
p3, = par2.plot([0, 1, 2], [50, 30, 15], label="Velocity")
host.set_xlim(0, 2)
host.set_ylim(0, 2)
par1.set_ylim(0, 4)
par2.set_ylim(1, 65)

host.set_xlabel("Distance")
host.set_ylabel("Density")
par1.set_ylabel("Temperature")
par2.set_ylabel("Velocity")


host.legend()

host.axis["left"].label.set_color(p1.get_color())
par1.axis["right"].label.set_color(p2.get_color())
par2.axis["right2"].label.set_color(p3.get_color())