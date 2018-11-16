# Strago
Interactive LC viewer and transit fit

![Strago](https://github.com/astrofelipe/Strago/raw/master/img/strago.png)
![Preview](https://github.com/astrofelipe/Strago/raw/master/img/stragopreview.png)

# Usage
From terminal:

bokeh serve --show rstrago.py --args <file with data>
  
Data should be in two columns: time and flux

## Optional args
All of them must go after the "--args" (from Bokeh)

--period <P> --t0 <t0>: This forces period and ephemeris if you know them previously, 
otherwise it will be taken from a BLS periodogram.

--sigmaclip: Does a very basic sigma clipping. Useful if your data has outliers, but transits
may be deleted too.
  
# Requirements
- Web browser
- Bokeh
- SciPy
- NumPy
- [Batman](https://github.com/lkreidberg/batman)
- [bls](https://github.com/dfm/python-bls)

# TODO
- Use different P and t0 for model and data. Initial idea was let them fixed so students can
play with the rest of parameters. Better precision on t0 when it's extracted from BLS.
- emcee + triangle panel above the periodogram
- Python 3 :(
