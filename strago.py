from bokeh.plotting import figure, show
from bokeh.layouts import gridplot, widgetbox, row, layout
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, TextInput, Toggle
from bokeh.io import curdoc
from scipy.ndimage.filters import median_filter
from astropy.stats import LombScargle
import bls
import numpy as np
import argparse
import batman

parser = argparse.ArgumentParser(description='Lightcurve tools')
parser.add_argument('File')
parser.add_argument('--period', type=float, default=None)
parser.add_argument('--t0', type=float, default=None)
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
#plot.line('t', 'trn', source=ndata, line_width=1, color='lime')

#BLS
pgram = figure(width=1000, height=200, x_range=[0,20])
blsre = bls.eebls(t, nf, np.ones(len(t)), np.ones(len(t)), 50000, 1/30., 1e-4, 250, 0.01, 0.15)
freqs = 1 / np.arange(1/30., 1/30. + 50000*1e-4, 1e-4)

blsda = ColumnDataSource(data=dict(per=freqs, pow=blsre[0]))
pgram.line('per', 'pow', source=blsda)

#GLS
glsf, glsp = LombScargle(t, f).autopower(minimum_frequency=1/20., maximum_frequency=1/0.1)
glsfig = figure(width=1000, height=200)
glsdat = ColumnDataSource(data=dict(per=1/glsf, pow=glsp))

glsfig.line('per', 'pow', source=glsdat)

#Phased lightcurve
pha = figure(width=400, height=200, x_range=[-.03,.03])
per = blsre[1] if args.period is None else args.period

inn = np.median([blsre[-2], blsre[-1]])
t0  = t[0] + per*(inn/250.) if args.t0 is None else args.t0

ph     = (t-t0) / per % 1.0
ph[ph>0.5] -= 1.0
offset = ph[np.argmin(nf)]
ph -= offset

phsrc = ColumnDataSource(data=dict(ph=ph, pf=nf))
pha.circle('ph', 'pf', source=phsrc, size=2, line_color='black')

#Batman
rprs   = 0.01
params = batman.TransitParams()
params.t0 = t0                       #time of inferior conjunction
params.per = per                      #orbital period
params.rp = np.sqrt(rprs)            #planet radius (in units of stellar radii)
params.a = 15.                       #semi-major axis (in units of stellar radii)
params.inc = 90.                     #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.1, 0.3]                #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model

tb = np.copy(t)#np.linspace(t.min(), t.max(), 10000)
pb = (tb - params.t0) / params.per % 1.0
pb[pb>0.5] -= 1.0
psort = np.argsort(pb)
mb = batman.TransitModel(params, tb)
fb = mb.light_curve(params)

bsrc = ColumnDataSource(data=dict(pb=pb[psort], fb=fb[psort], tb=tb, fb2=fb*trn))
pha.line('pb', 'fb', source=bsrc, line_color='firebrick')

#Widgets BATMAN
brp  = Slider(title='Radio planeta (R_Earth)', value=10, start=0.5, end=50, step=1e-5)
brs  = Slider(title='Radio estrella (R_Sun)', value=1, start=0.05, end=10, step=1e-5)
ba   = Slider(title='Distancia a la estrella (UA)', value=0.05, start=0.001, end=1.5, step=1e-4)
becc = Slider(title='Excentricidad', value=params.ecc, start=0, end=1, step=1e-3)
binc = Slider(title='Inclinación (grados)', value=params.inc, start=80, end=100, step=1e-5)
bu1  = Slider(title='Limb Darkening u1', value=params.u[0], start=0, end=1, step=1e-4)
bu2  = Slider(title='Limb Darkening u2', value=params.u[1], start=0, end=1, step=1e-4)

#Widgets "basic"
#period = Slider(title='Period', value=per, start=0.1, end=40., step=0.001)
#t0w    = Slider(title='t0', value=t0, start=t0-2, end=t0+2, step=1e-6)
b = params.a * np.cos(np.radians(params.inc))

txtper = TextInput(value=str(per), title='Periodo')
txtt0  = TextInput(value=str(t0), title='t0')
txtb   = TextInput(value='%.3f' % b, title='Parametro de impacto')


inputs = widgetbox(txtper, txtt0, txtb)

def update_batman(attrname, old, new):
    rps = brp.value * 0.009158 #Earth radii in Sun radii
    ars = ba.value * 215 #AU to Sun radii

    #params.per = bper.value
    params.rp  = rps/brs.value
    #params.t0  = bt0.value
    params.a   = ars*brs.value
    params.u   = [bu1.value, bu2.value]
    params.ecc = becc.value
    params.inc = binc.value

    #pb = ((tb - t0.value) / bper.value) % 1.0
    #pb[pb>0.5] -= 1.0
    #psort = np.argsort(pb)
    fb = mb.light_curve(params)
    bsrc.data = dict(pb=pb[psort], fb=fb[psort], tb=tb, fb2=fb*trn)

    #resi = (fb - nf)
    #resrc.data = dict(ph=pb, resi=resi)

    er = (ba.value*(1-becc.value**2))/(1+becc.value*np.cos(et))
    ex = er*np.cos(et)
    ey = er*np.sin(et)
    odata.data = dict(x=ex, y=ey)
    sundat.data = dict(xc=[0], yc=[0], s=[brs.value], sn=[brs.value/215.])

    b = params.a * np.cos(np.radians(params.inc))
    txtb.value  = '%.3f' % b
    pladat.data = dict(xc=[0], yc=[b], s=[brp.value*0.009158])

    '''
    xor = np.linspace(-2.5*brs.value, 2.5*brs.value, 100)
    yor = np.zeros(100)

    angl   = np.radians(params.inc) - np.pi/2.
    rotmat = np.array([[np.cos(angl), np.sin(angl)],[-np.sin(angl), np.cos(angl)]])

    xor, yor = np.dot(rotmat, np.array([xor, yor]))

    orbdat.data = dict(x=xor, y=yor)
    '''

for w in [brp, brs, ba, becc, binc, bu1, bu2]:
    w.on_change('value', update_batman)

binput = widgetbox(brp, brs, ba, becc, binc, bu1, bu2)

'''
#Residuals
phar = figure(width=400, height=75, x_range=[-.03,.03])
resi = (fb - nf)
resrc = ColumnDataSource(data=dict(ph=pb, resi=resi))
phar.circle('ph', 'resi', source=resrc, size=2, line_color='black')
'''

#Orbit plot
orbit = figure(width=400, height=400, match_aspect=True)

et = np.linspace(0, 2*np.pi, 360)
er = (ba.value*(1-becc.value**2))/(1+becc.value*np.cos(et))
ex = er*np.cos(et)
ey = er*np.sin(et)

odata = ColumnDataSource(data=dict(x=ex, y=ey))
orbit.line('x', 'y', source=odata)

sundat = ColumnDataSource(data=dict(xc=[0], yc=[0], s=[brs.value], sn=[brs.value/215.]))
osun   = orbit.circle('xc', 'yc', radius='sn', source=sundat, fill_color='yellow')


#Line of sight view
los = figure(width=400, height=400, x_range=[-3,3], y_range=[-3,3])
#orbdat = ColumnDataSource(data=dict(x=np.linspace(-2.5*brs.value, 2.5*brs.value, 100), y=np.zeros(100)))

los.circle('xc', 'yc', radius='s', source=sundat, fill_color='yellow')
#los.line('x','y', source=orbdat)

pladat = ColumnDataSource(data=dict(xc=[0], yc=[b], s=[brp.value*0.009158]))
los.circle('xc', 'yc', radius='s', source=pladat, fill_color='brown')




#All together now (8)
#curdoc().add_root(row(plot, pha, inputs))
#gplot = gridplot([[plot, None],[inputs, binput, pha],[None, None, phar],[pgram, None],[glsfig,None]])
laylay = layout([[plot, pha],[inputs, binput, orbit, los],[pgram]])
curdoc().add_root(laylay)
curdoc().title = 'Lightcurve'

#pp = gridplot([[plot, pha, inputs]])
#show(pp)
