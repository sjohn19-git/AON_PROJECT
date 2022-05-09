import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize

import obspy
from obspy.core.util import AttribDict
from obspy.imaging.cm import obspy_sequential
from obspy.signal.invsim import corn_freq_2_paz
from obspy.signal.array_analysis import array_processing

import matplotlib.dates as mdates
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

t1 = UTCDateTime("2019-11-13T01:01:00.000")
t2 = UTCDateTime("2019-11-13T01:02:00.000")
fdsn_client = Client('IRIS')
# Fetch waveforms and reponse from IRIS FDSN web service into ObsPy
# stream and inventory objects

st = fdsn_client.get_waveforms(
    network='TA', station='M16K', location='*', channel='BHZ',
    starttime=t1, endtime=t2)
inv = fdsn_client.get_stations(
    network='TA', station='M16K', location='*', channel='BHZ',
    starttime=t1, endtime=t2, level='response')
# define a filter band to prevent amplifying noise during the deconvolution
pre_filt = (0.08,0.1,0.2,0.22)
st.remove_response(inventory=inv, output='DISP', pre_filt=pre_filt,taper=True)
st.detrend(type="demean")


st0=st #add trace to new stream

st = fdsn_client.get_waveforms(
    network='TA', station='M17K', location='*', channel='BHZ',
    starttime=t1, endtime=t2)
inv = fdsn_client.get_stations(
    network='TA', station='M17K', location='*', channel='BHZ',
    starttime=t1, endtime=t2, level='response')
# define a filter band to prevent amplifying noise during the deconvolution
pre_filt = (0.005, 0.006, 30.0, 35.0)
st.remove_response(inventory=inv, output='DISP', pre_filt=pre_filt)
st.detrend(type="demean")
st0+=st

st = fdsn_client.get_waveforms(
    network='TA', station='N16K', location='*', channel='BHZ',
    starttime=t1, endtime=t2)
inv = fdsn_client.get_stations(
    network='TA', station='N16K', location='*', channel='BHZ',
    starttime=t1, endtime=t2, level='response')
# define a filter band to prevent amplifying noise during the deconvolution
pre_filt = (0.08,0.1,0.2,0.22)
st.remove_response(inventory=inv, output='DISP', pre_filt=pre_filt)
st.detrend(type="demean")
st0+=st

st = fdsn_client.get_waveforms(
    network='TA', station='N17K', location='*', channel='BHZ',
    starttime=t1, endtime=t2)
inv = fdsn_client.get_stations(
    network='TA', station='N17K', location='*', channel='BHZ',
    starttime=t1, endtime=t2, level='response')
# define a filter band to prevent amplifying noise during the deconvolution
pre_filt = (0.08,0.1,0.2,0.22)
st.remove_response(inventory=inv, output='DISP', pre_filt=pre_filt)
st.detrend(type="demean")
st0+=st

st0.plot()




st0[0].stats.coordinates = AttribDict({
    'latitude': 61.0224,
    'elevation': 0.39,
    'longitude': 180+180-158.9593})

st0[1].stats.coordinates = AttribDict({
    'latitude': 61.4009,
    'elevation': 0.387,
    'longitude': 180+180-157.4375})

st0[2].stats.coordinates = AttribDict({
    'latitude': 60.4742,
    'elevation': 0.582,
    'longitude': 360-158.769})

st0[3].stats.coordinates = AttribDict({
    'latitude': 60.5269,
    'elevation': 0.350,
    'longitude': 360-157.1867})









stime = t1+1
etime = t2-1

kwargs = dict(
    # slowness grid: X min, X max, Y min, Y max, Slow Step
    sll_x=-3, slm_x=3, sll_y=-3, slm_y=3, sl_s=0.03,
    # sliding window properties
    win_len=1.0, win_frac=0.05,
    # frequency properties
    frqlow=0.2, frqhigh=0.4, prewhiten=0,
    # restrict output
    semb_thres=-1e9, vel_thres=-1e9, timestamp='mlabday',
    stime=stime, etime=etime
)
out = array_processing(st0, **kwargs)


#plotting--------------------------------------

labels = ['rel.power', 'abs.power', 'baz', 'slow']

xlocator = mdates.AutoDateLocator()
fig = plt.figure(figsize=(6,6))
for i, lab in enumerate(labels):
    ax = fig.add_subplot(4, 1, i + 1)
    ax.scatter(out[:, 0], out[:, i + 1], c=out[:, 1], alpha=0.6,
               edgecolors='none', cmap=obspy_sequential)
    ax.set_ylabel(lab)
    ax.set_xlim(out[0, 0], out[-1, 0])
    ax.set_ylim(out[:, i + 1].min(), out[:, i + 1].max())
    ax.xaxis.set_major_locator(xlocator)
    ax.xaxis.set_major_formatter(mdates.AutoDateFormatter(xlocator))

fig.suptitle('Microseism'+stime.strftime('%Y-%m-%d'))
fig.autofmt_xdate()
fig.subplots_adjust(left=0.15, top=0.95, right=0.95, bottom=0.2, hspace=0)
plt.show()






cmap = obspy_sequential

# make output human readable, adjust backazimuth to values between 0 and 360
t, rel_power, abs_power, baz, slow = out.T
baz[baz < 0.0] += 360

# choose number of fractions in plot (desirably 360 degree/N is an integer!)
N = 36
N2 = 30
abins = np.arange(N + 1) * 360. / N
sbins = np.linspace(0, 3, N2 + 1)

# sum rel power in bins given by abins and sbins
hist, baz_edges, sl_edges = \
    np.histogram2d(baz, slow, bins=[abins, sbins], weights=rel_power)

# transform to radian
baz_edges = np.radians(baz_edges)

# add polar and colorbar axes
fig = plt.figure(figsize=(8, 8))
cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
ax.set_theta_direction(-1)
ax.set_theta_zero_location("N")

dh = abs(sl_edges[1] - sl_edges[0])
dw = abs(baz_edges[1] - baz_edges[0])

# circle through backazimuth
for i, row in enumerate(hist):
    bars = ax.bar((i * dw) * np.ones(N2),
                  height=dh * np.ones(N2),
                  width=dw, bottom=dh * np.arange(N2),
                  color=cmap(row / hist.max()))

ax.set_xticks(np.linspace(0, 2 * np.pi, 4, endpoint=False))
ax.set_xticklabels(['N', 'E', 'S', 'W'])

# set slowness limits
ax.set_ylim(0,3)
[i.set_color('grey') for i in ax.get_yticklabels()]
ColorbarBase(cax, cmap=cmap,
             norm=Normalize(vmin=hist.min(), vmax=hist.max()))

plt.show()


