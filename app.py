import io
import streamlit as st
import matplotlib.pyplot as plt
from pycbc.waveform import get_td_waveform
from pycbc.detector import Detector
import numpy as np
from scipy import signal
from scipy.io import wavfile
from gwpy.timeseries import TimeSeries

np.random.seed(38)

def make_audio_file(bp_data, t0=None):
    # -- window data for gentle on/off
    window = signal.windows.tukey(len(bp_data), alpha=1.0/100)
    win_data = bp_data*window

    # -- Normalize for 16 bit audio
    win_data = np.int16(win_data/np.max(np.abs(win_data)) * 32767 * 0.9)

    fs=1/win_data.dt.value
    virtualfile = io.BytesIO()    
    wavfile.write(virtualfile, int(fs), win_data)
    
    return virtualfile


st.set_page_config(page_title="Waveform analysis", page_icon=":scroll:", layout="centered")
st.title("Waveform analysis")

recipelist = ["GW150914", "GW170817"]

st.sidebar.header("🪄 Waveform wizard")

chosen = st.sidebar.selectbox("Select Recipe", recipelist)
st.subheader(chosen.upper())

    

if chosen=="GW150914":
    apx = 'SEOBNRv4'
    ra, dec = 2.7, -1.2
    f_lower = 20
    incl = 140.
    m1D, m2D = 35. , 30. 
    dLD = 470.
    iD  = 2.45
    phaseD = 0.5
    xlim = [-0.4,0.05]

if chosen=="GW170817":
    apx = 'IMRPhenomD'
    ra, dec = 3.3, 0.3
    f_lower = 50
    incl = 100.
    m1D, m2D = 1.45 , 1.28 
    dLD = 40.
    iD  = 2.45
    phaseD = 0.9
    xlim = (None, 0)
    xlim = [-1.,0.05]

  



col1, col2 = st.sidebar.columns(2)

with col1:
    m1 = col1.slider('$m_1~\mathrm{(M_\odot)}$', 1., 100.0, 10.0, step=0.01)
with col2:
    m1 = col2.number_input(" ", value=m1)

col1, col2 = st.sidebar.columns(2)
with col1:
    m2 = col1.slider('$m_2~\mathrm{(M_\odot)}$', 1., 100.0, 10.0, step=0.01)
with col2:
    m2 = col2.number_input("  ", value=m2)

col1, col2 = st.sidebar.columns(2)
with col1:
    dL = col1.slider('$d_L~\mathrm{(Mpc)}$', 0.01, 500., 50., step=1.)
with col2:
    dL = col2.number_input("   ", value=dL)

col1, col2 = st.sidebar.columns(2)
with col1:
    phase = col1.slider('$\phi~\mathrm{(deg)}$', 0., 3.14, 0.9, step=0.1)
with col2:
    phase = col2.number_input("    ", value=phase)

col1, col2 = st.sidebar.columns(2)
with col1:
    incl  = col1.slider('$\iota~\mathrm{(deg)}$', 0., 360., incl, step=1.) 
with col2:
    incl = col2.number_input("     ", value=incl)



xlim[0] = st.sidebar.slider('$x_{min}$', -30., 0., xlim[0], step=.1) 


# NOTE: Inclination runs from 0 to pi, with poles at 0 and pi
#       coa_phase runs from 0 to 2 pi.
hp, hc = get_td_waveform(approximant=apx,
                         mass1=m1,
                         mass2=m2,
                         distance=dL,
                         inclination=incl / 180 * np.pi,
                         coa_phase=phase,
                         delta_t=1.0/4096,
                         f_lower=f_lower
                         )

det_l1 = Detector('L1')


hpData, hcData = get_td_waveform(approximant=apx,
                                 mass1=m1D,
                                 mass2=m2D,
                                 distance=dLD,
                                 inclination= iD / 180 * np.pi,
                                 coa_phase=phaseD,
                                 delta_t=1.0/4096,
                                 f_lower=f_lower
                                 )


end_time = 0
declination = dec
right_ascension = ra
polarization = 2.34
hp.start_time += end_time
hc.start_time += end_time

signal_l1     = det_l1.project_wave(hp, hc,  right_ascension, declination, polarization)
signal_l1Data = det_l1.project_wave(hpData, hcData,  right_ascension, declination, polarization)


fig = plt.figure(figsize=(7,4), dpi=200)

y = signal_l1Data 
y = y + 0.5*np.random.normal(loc=0, scale=np.std(y), size=len(y))

plt.plot(signal_l1Data.sample_times, y, label='LIGO Livingston', lw=1.5, c='cornflowerblue')
plt.plot(signal_l1.sample_times, signal_l1, label='Model', lw=.8, c='k')

plt.ylabel('Strain')
plt.xlabel('Time (s)')
plt.legend()
plt.xlim(*xlim)

st.pyplot(plt)


ts = TimeSeries(y, dt=1.0/4096).taper() 
tm = TimeSeries(signal_l1, dt=1.0/4096).taper() 
st.write("Sonification of data")
st.audio(make_audio_file(ts), format='audio/wav')
st.write("Sonification of model")
st.audio(make_audio_file(tm), format='audio/wav')

table = [[r"$m_1$", "Mass of the primry object"]]



st.markdown("""
| Parameter | Description |
|-----------|----------------|
| $m_1$     | Mass of the primary object (in solar masses) |
| $m_2$     | Mass of the secondary object (in solar masses) |
| $d_L$     | Distance to the binary (in Mpc) |
| $\phi$    | Coalesence phase of the binary (in degrees) |
| $\iota$   | Angle between orbital angular momentum L and line-of-sight at the reference frequency |
""")