import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from WF4Py import waveforms as WF

dataET    = np.genfromtxt("sensitivity_ET.dat")
dataLVK   = np.genfromtxt("sensitivity_LVK.dat")
dataaGW15 = np.genfromtxt("GW150914.txt")

def Mc(m1,m2):    
    return (m1*m2)**(3/5) / (m1+m2)**(1/5)

def eta(m1,m2):
    return (m1*m2) / (m1+m2)**2

st.title('Plot Model (IMRPhenomD)')

m1 = 10
m2 = 10
dL = 0.01

tmpWF = WF.IMRPhenomD()





sectionnames = [
                'Time Domain',
                'Frequency Domain',
]

def headerlabel(number):
    return "{0}: {1}".format(number, sectionnames[number-1])
    
page = st.radio('Sezione:', [1,2], format_func=headerlabel)

if page == 1:

    m1 = st.slider('Mass 1', 1.2, 100.0, 10.0, step=0.1)
    m2 = st.slider('Mass 2', 1.2, 100.0, 10.0, step=0.1)
    dL = st.slider('Distance in Mpc', 10, 1000, 40)

    events = {'Mc':np.array([Mc(m1,m2)]),
              'dL':np.array([dL/1000]),
              'iota':np.array([0.]),
              'eta':np.array([0.24]),
              'chi1z':np.array([0.8]),
              'chi2z':np.array([-0.8]),
              'Lambda1':np.array([0.]),
              'Lambda2':np.array([0.])}

    fcut = tmpWF.fcut(**events)
    fminarr = np.full(fcut.shape, 1)
    fgrids = np.geomspace(fminarr, fcut, num=int(1000))

    myampl = tmpWF.Ampl(fgrids, **events)
    myphase = tmpWF.Phi(fgrids, **events)
    mytime = tmpWF.tau_star(fgrids, **events)

    fig_time = plt.figure(figsize=(8,5))
    # Number of samples in the time-series representation

    # Compute the time-series representation of the signal using IDFT
    N = 10000
    x = np.fft.ifft(myampl.flatten() + 1j * fgrids.flatten(), N)
    t = np.arange(0,1,1/N)


    plt.plot(dataaGW15[:,0],dataaGW15[:,1])
    # plt.plot(10*(t[9500:]-1),x[9500:]*1e-22, ls="-", c='k', label="Model")
    # plt.loglog()
    plt.ylabel('Amplitude')
    plt.xlabel('Time (s)')
    # # plt.axis('off')
    # plt.tick_params(left = False, right = False , labelleft = False ,
    #                 labelbottom = False, bottom = False)

    st.markdown("")
    st.markdown("""**Waveform** (Amplitude vs. Time) with arbitrary units""")
    st.pyplot(fig_time)


          
if page == 2:

    m1 = st.slider('Mass 1', 1.2, 100.0, 10.0, step=0.1)
    m2 = st.slider('Mass 2', 1.2, 100.0, 10.0, step=0.1)
    dL = st.slider('Distance in Mpc', 10, 2000, 40)
    mc = Mc(m1,m2)*1.1
    events = {'Mc':np.array([Mc(m1,m2)*1.1]),
            'dL':np.array([dL/1000]),
            'iota':np.array([0.]),
            'eta':np.array([0.24]),
            'chi1z':np.array([0.8]),
            'chi2z':np.array([-0.8]),
            'Lambda1':np.array([0.]),
            'Lambda2':np.array([0.])}

    fcut = tmpWF.fcut(**events)
    fminarr = np.full(fcut.shape, 1)
    fgrids = np.geomspace(fminarr, fcut, num=int(1000))

    myampl = tmpWF.Ampl(fgrids, **events)
    myphase = tmpWF.Phi(fgrids, **events)
    mytime = tmpWF.tau_star(fgrids, **events)

    st.markdown("Chirp mass: "+str(round(mc,2)))

    fig = plt.figure(figsize=(8,5))
    plt.plot(fgrids, myampl, ls="--", color='k', label="Model")
    plt.plot(dataLVK[:,0],dataLVK[:,1],color='teal', label="aLIGO")
    plt.plot(dataET[:,0],dataET[:,3],color='goldenrod', label="ET")
    plt.grid(True, which='both', ls='dotted', linewidth='0.8', alpha=.8)
    plt.ylabel('Strain')
    plt.xlabel('Frequency (Hz)')
    plt.ylim(5e-26,1e-17)
    plt.legend()
    plt.loglog()
    st.pyplot(fig)











