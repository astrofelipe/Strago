from bokeh.plotting import figure, show
from bokeh.layouts import gridplot, widgetbox, row
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput, Toggle
from bokeh.io import curdoc
from scipy.ndimage.filters import median_filter
from astropy.stats import LombScargle
import bls
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Lightcurve tools')
parser.add_argument('File')
parser.add_argument('--period', type=float, default=1.0)
parser.add_argument('--t0', type=float, default=0.0)
args = parser.parse_args()

#Data
fname = args.File
t,f   = np.genfromtxt(fname, unpack=True, usecols=(0,1))
src   = ColumnDataSource(data=dict(t=t, f=f))

#Normalize
Nmf = int(np.sqrt(len(t)))
trn = median_filter(f, size=Nmf)
nf  = f / trn
ndata = ColumnDataSource(data=dict(t=t, trn=trn))

#Full lightcurve
plot = figure(plot_height=200, plot_width=1000, title='Lightcurve',
              x_range=[np.nanmin(t), np.nanmax(t)])

plot.circle('t', 'f', source=src, size=1)
plot.line('t', 'trn', source=ndata, line_width=1, color='lime')

#BLS
pgram = figure(width=1000, height=200, x_range=[0,20])
blsre = bls.eebls(t, nf, np.ones(len(t)), np.ones(len(t)), 50000, 1/30., 1e-4, 200, 0.01, 0.15)
freqs = 1 / np.arange(1/30., 1/30. + 50000*1e-4, 1e-4)

blsda = ColumnDataSource(data=dict(per=freqs, pow=blsre[0]))
pgram.line('per', 'pow', source=blsda)

#GLS
glsf, glsp = LombScargle(t, f).autopower(minimum_frequency=1/20., maximum_frequency=1/0.1)
glsfig = figure(width=1000, height=200)
glsdat = ColumnDataSource(data=dict(per=1/glsf, pow=glsp))

glsfig.line('per', 'pow', source=glsdat)

#Phased lightcurve
pha = figure(width=400, height=200, x_range=[-.05,.05])
per = blsre[1]#args.period

inn = np.mean([blsre[-2], blsre[-1]])
t0  = t[0] + per*(inn/200.)#args.t0

ph    = (t-t0) / per % 1.0
ph[ph>0.5] -= 1.0
phsrc = ColumnDataSource(data=dict(ph=ph, pf=nf))
pha.circle('ph', 'pf', source=phsrc, size=1)


#Widgets
period = Slider(title='Period', value=per, start=0.1, end=40., step=0.001)
t0w    = Slider(title='t0', value=t0, start=0.0, end=4000, step=0.001)
txtper = TextInput(value=str(per), title='Period')


#Callbacks
def update_data(attrname, old, new):
    pf  = phsrc.data['pf']
    per = period.value
    t0  = t0w.value

    ph    = (t-t0) / per % 1.0
    ph[ph>0.5] -= 1.0

    phsrc.data = dict(ph=ph, pf=pf)

for w in [period, t0w]:
    w.on_change('value', update_data)

inputs = widgetbox(period, t0w, txtper)

#All together now (8)
#curdoc().add_root(row(plot, pha, inputs))
curdoc().add_root(gridplot([[plot, pha],[pgram, None],[glsfig,None],[None, inputs]]))
curdoc().title = 'Lightcurve'

#pp = gridplot([[plot, pha, inputs]])
#show(pp)
