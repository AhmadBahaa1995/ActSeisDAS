# Those are libraries and function i usually call, it doesn't haram if you call it and not use it

import mpu
import obspy
import sys
import csv
import scipy
import os
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import time
import math
import multiprocessing

from obspy import UTCDateTime, read, Stream, Trace
from obspy.signal import cross_correlation
from obspy.signal.polarization import polarization_analysis
from obspy.signal.util import next_pow_2

from numpy.fft import fft,ifft,fftfreq
from numpy import var

from scipy.interpolate import interp1d
from scipy.fftpack import rfft, irfft, fftfreq
from scipy.signal import butter, lfilter, correlate, coherence, minimum_phase, hilbert, find_peaks
from scipy.stats import spearmanr, pearsonr, kendalltau
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import shift

from datetime import date, datetime
from IPython.display import clear_output
from time import time
from ast import While
from multiprocessing import Pool

def t_nmo(offset, t0=0.2, v_nmo=4000):
    '''
    Calculate NMO Time
    
    Input:
    offset = Distance location of each trace [m]. Type: 1-D list or array.
    t0 = Time at offset=0 or starting time of the trace [s]. Type: float. Default: 0.2.
    v_nmo = NMO velocity that shaped the curve [m/s]. Type: float. Default: 4000.
    
    Output:
    Arrival time of the trace [s]. Type: 1-D array with same size as offset.
    '''
    return np.sqrt(t0**2 + (offset**2)/(v_nmo**2))
	
def RickerWavelet(f=25, l=0.128, dt=0.002):
    '''
    Create Ricker (Mexican hat) wavelet to be convoluted with the reflectivity.
    
    Input:
    f = Wavelet frequency [Hz]. Lower is wider. Type: float. Default: 25.
    l = Time length [s]. Type: float. Default: 0.128.
    dt = Sampling rate [s]. Type: float. Default: 0.002 (2 ms).
    
    Output:
    t = Wavelet time [s]. Type: 1-D array.
    y = Wavelet amplitude. Maximum 1. Minimum < 0. Type: 1-D array.
    '''
    t = np.arange(-l/2, (l-dt)/2, dt)
    y = (1 - 2*((np.pi)**2)*(f**2)*(t**2)) / np.exp(((np.pi)**2)*(f**2)*(t**2))
    return t, y

def create_synthetic(time, velocity, channel, offset, dt=0.002, wavelet_freq=25):
    '''
    Create synthetic seismogram from NMO time. Convert the time to the reflectivity amplitude,
    then, convoluted with the wavelet.
    
    Input:
    time = Zero-offset time of the trace [s]. Type: 1-D list or array.
    velocity = NMO velocity of the curve [m/s]. Type: 1-D list or array.
    ***Length of time and velocity are MUST the same. It represents the layers.***
    channel = Number of channel (trace). Type: int.
    dt = Sampling rate [s]. Type: float. Default: 0.002 (2 ms).
    wavelet_freq = Wavelet frequency [Hz]. Lower is wider. Type: float. Default: 25.
    
    Output:
    amplitude = Amplitude matrix of each traces over time. Type: 2-D array with size (time/dt) x channel.
    Y = Time coordinate of one trace [s]. Type: 1-D array with size time/dt.
    '''
    t = np.zeros((len(time), channel))
    vel = np.zeros((len(velocity), channel))
    
    for i in range(len(vel)):
        vel[i].fill(velocity[i])
        t[i].fill(time[i])
        
    NMO = t_nmo(offset, t0=t, v_nmo=vel)
    
    X = np.arange(0,channel,1)
    Y = np.arange(0, np.max(NMO).round(), dt)
    time_matrix = np.meshgrid(X,Y)[1]
    
    target = np.array([np.abs(Y-i).argmin() for i in NMO.flatten()]).reshape(NMO.shape)
    
    amplitude_matrix = np.zeros(time_matrix.shape)
    
    for i in range(NMO.shape[1]):
        time_matrix[target[:,i],i] = NMO[:,i]
        amplitude_matrix[target[:,i],i] = 1
        
    t1, y1 = RickerWavelet(f=wavelet_freq)
    
    amplitude = np.array([np.convolve(amplitude_matrix[:,i],y1,'same') for i in range(amplitude_matrix.shape[1])]).T
    
    return amplitude, Y

def plot_seismogram(ax, data, offset, time, plot_option=1, cmap='seismic'):
    '''
    Simplified function to easily plot the seismogram.
    
    Input:
    ax = matplot.pyplot.Axes of the plot. Might be declared before or not, depends on the plot_option. 
    data = Seismogram amplitude matrix data over the time. Type: 2-D array.
    offset = Distance location of each trace [m]. Type: 1-D list or array with size of data.shape[1].
    time = Time coordinate of one trace [s]. Type: 1-D array with size of data.shape[0].
    plot_option = Plot view option of the seismogram. Explained in the function comment below. Type: int or str. Default: 1.
    cmap = Colormap of the seismogram. Not necessary if plot_option = 0 or 3 or 'wiggle' or 'stacked'. Type: str.
    ***Further options of cmap : https://matplotlib.org/stable/tutorials/colors/colormaps.html***
    
    If the plot_option is not 2 or 'both', the ax must be declared before.
    Example:
    fig, ax = plt.subplots()
    plot_seismogram(ax, your_data, your_offset, your_time, plot_option=1, cmap='seismic')
    
    Output:
    matplotlib.pyplot show of the seismogram plot (no return).
    '''
    ### Wiggle trace
    def wiggle(ax, data, time, plot_option=plot_option):
        for i in range(amplitude.shape[1]):
            ax.plot(amplitude[:,i]+i, Y, linewidth=0.2, c='k')
            ax.fill_between(amplitude[:,i]+i,Y,0, where=(amplitude[:,i]+i>=i) ,lw=0, color='k')
        
        if plot_option == 0 or plot_option == 'wiggle':
            ax.invert_yaxis()
        
        ax.set_ylabel('Time [s]')
        ax.set_xlabel('Trace number')
        
    ### Filled color trace
    def color(ax, data, time, offset, plot_option=plot_option):
        ax.imshow(data, cmap=cmap,
                  extent=[np.min(offset), np.max(offset),
                          np.max(time), np.min(time)],
                  aspect='auto', vmin=-100, vmax=100)
        
        if plot_option == 1 or plot_option == 'color':
            ax.set_ylabel('Time [s]')
            
        ax.set_xlabel('Offset [m]')
    
    ### Stacked (single) trace of the data input
    def stack_seismogram(ax, data, time, offset, plot_option=plot_option):
        stacked = data.sum(axis=1)
        
        positive_amp = stacked.copy()
        positive_amp[positive_amp<0] = 0
        
        ax.plot(stacked, time, lw=0.5, color='k')
        ax.fill_between(positive_amp, time, color='k', lw=0)
        ax.set_xlabel('Amplitude')
    
    # Plot the wiggle trace only. ax must be declared before
    if plot_option == 0 or plot_option == 'wiggle':
        wiggle(ax, data, time)
    # Plot the filled color trace only. ax must be declared before
    elif plot_option == 1 or plot_option == 'color':
        color(ax, data, time, offset)
    # Plot both of wiggle and filled color trace. ax does not need to be declared before.
    elif plot_option == 2 or plot_option == 'both':
        fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True)
        wiggle(ax[0], data, time)
        color(ax[1], data, time, offset)
    # Plot the stacked trace only. ax must be declared before
    elif plot_option == 3 or plot_option == 'stacked':
        stack_seismogram(ax, data, time, offset)

def shift_nmo(data, t0, v, time, offset):
    '''
    Shifting the seismogram matrix based on the velocity and time from NMO.
    This code was modified from "Geophysics I: Theory of Geophysical Prospection Methods" exercise course
    by Prof. Dr.sc. Florian M. Wagner, Geophysical Imaging and Monitoring (GIM) RWTH Aachen University.
    
    Input:
    data = Seismogram amplitude matrix data over the time. Type: 2-D array.
    t0 = Time at offset=0 or starting time of the trace [s]. Type: float.
    v = Velocity of the curve [m/s]. Type: float.
    time = Time coordinate of one trace [s]. Type: 1-D array with the size of data.shape[0].
    offset = Distance location of each trace [m]. Type: 1-D list or array with the size of data.shape[1].
    
    Output:
    zeros_nmo = Shifted data based on v and t0. Type: 2-D array with the same size of data.
    '''
    zeros_nmo = np.zeros_like(data)
    i = 0
    for trace in data.T:
        nmo = t_nmo(offset[i], t0=t0, v_nmo=v)
        new_nmo = nmo - (t0)
        target_nmo = np.abs(new_nmo - time).argmin()
        zeros_nmo[:,i] = shift(trace, -target_nmo)
        i += 1
    return zeros_nmo

def NMO_correction(data,offset,time,v=4000, t0=0.2):
    '''
    Visualization of the NMO correction based on the shifting.
    This code was modified from "Geophysics I: Theory of Geophysical Prospection Methods" exercise course
    by Prof. Dr.sc. Florian M. Wagner, Geophysical Imaging and Monitoring (GIM) RWTH Aachen University.
    
    This is interactive visualization to play around with the zero-offset time (t0) and the correct
    NMO velocity (v)
    
    Input:
    data = Seismogram amplitude matrix data over the time. Type: 2-D array.
    offset = Distance location of each trace [m]. Type: 1-D list or array with size of data.shape[1].
    time = Time coordinate of one trace [s]. Type: 1-D array with size of data.shape[0].
    t0 = Time at offset=0 or starting time of the trace [s]. Type: float. Default: 0.2.
    v = Velocity of the curve [m/s]. Type: float. Default: 4000.
    
    Output:
    Four matplotlib.pyplot plots with:
    Top left = data with the NMO time plot (to be matched with the curve).
    Top right = Stacked data without applying shift. Always the same.
    Bottom left = Shifted data based on choosen v and t0.
    Bottom right = Stacked data after shifting.
    '''
    fig, ax = plt.subplots(nrows=2,ncols=2, figsize=(14,12), gridspec_kw={'width_ratios': [3, 1]}, sharey=True)        
    
    def velocity_test(ax,data,offset,time,v=v, t0=t0):
        nmo = t_nmo(offset, t0=t0, v_nmo=v)
        plot_seismogram(ax, data, offset, time, plot_option=1)
        ax.plot(offset, nmo, c='yellow', linewidth=5)
        
    velocity_test(ax[0,0], data, offset, time)
    
    shifted_nmo = shift_nmo(data, t0, v, time, offset)
        
    plot_seismogram(ax[1,0], shifted_nmo, offset, time, plot_option=1)
    
    ax[0,1].set_title('Trace stacking',fontsize=15,fontweight='bold')
    
    plot_seismogram(ax[0,1], data, offset, time, plot_option=3)
    plot_seismogram(ax[1,1], shifted_nmo, offset, time, plot_option=3)
    
    ax[0,0].set_title('Before NMO',fontsize=15,fontweight='bold')
    ax[1,0].set_title('After NMO',fontsize=15,fontweight='bold')
	
def semblance(data, time, v_range, offset):
    '''
    Calculate semblance as a guidance for velocity analysis. Using formula based on Geldart & Sheriff (2004).
    The highest value will likely shows the correct velocity at certain time.
    
    ***Warning: Running so slow and so far still bad result when tested with the real data***
    
    Input:
    data = Seismogram amplitude matrix data over the time. Type: 2-D array.
    time = Time coordinate of one trace [s]. Type: 1-D array with size of data.shape[0].
    v_range = Range of velocity to be scanned each [m/s]. More is slower, but smoother result. Range from lowest to highest. Type: 1-D list or array.
    offset = Distance location of each trace [m]. Type: 1-D list or array with size of data.shape[1].
    
    Output:
    Semblance = Coherency matrix of the summed trace amplitude based on velocity over the time. Type: 2-D array with the size of time x v_range.
    '''
    Semblance = np.zeros(( len(time), len(v_range) ))
    i = 0
    for T0 in time:
        j = 0
        for vel in v_range:
            semb = shift_nmo(data, round(T0, 2), vel, time, offset)
            Semblance[i,j] = ( (semb[i,:].sum())**2 ) / ( len(semb[i,:])*(semb[i,:]**2).sum() )
            j += 1
        i += 1
    return Semblance

def semblance_analysis(Semblance, data, t0, time, offset, vel_range, v=2000, Veloo=2000, Timee=0.2):
    '''
    Visualization of the semblance result, data seismogram, and the stacked. Interactive button
    and slider should be from different function because so far it still hard to put it here.
    
    Input:
    Semblance = semblance matrix result. Type: 2-D array.
    data = Seismogram amplitude matrix data over the time. Type: 2-D array.
    t0 = widgets.FloatSlider of the zero-offset time [s]. Used in the interactive mode. Type: float.
    time = Time coordinate of one trace [s]. Type: 1-D array with size of data.shape[0].
    offset = Distance location of each trace [m]. Type: 1-D list or array with size of data.shape[1].
    v = widgets.FloatSlider of the NMO velocity [m/s]. Used in the interactive mode. Type: float. Default: 2000.
    Veloo = Should be empty list ([]) in the interactive mode to store the NMO velocity from v. But here declared the default because
            sometimes it's error without it.
    Timee = Should be empty list ([]) in the interactive mode to store the zero-offset time from t0. But here declared the default because
            sometimes it's error without it.
    
    Output:
    Three matplotlib.pyplot plots with:
    Left = Semblance panel to interact and choose the best velocity over the time.
    Center = Seismogram amplitude matrix to be shifted based on the semblance picked.
    Right = Stacked seismogram after shifted based on the semblance picked.
    '''
    Shift_NMO = shift_nmo(data, t0, v, time, offset)
    
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(14,10), gridspec_kw={'width_ratios': [2, 3, 1]}, sharey=True)

    ax[0].imshow(Semblance,aspect='auto', cmap='jet',
                   extent=[np.min(vel_range), np.max(vel_range),
                           np.max(time), np.min(time)])
    ax[0].set_title('Semblance', fontsize=14, fontweight='bold')
    ax[0].set_ylabel('Time [s]')
    ax[0].set_xlabel('Velocity [m/s]')
    ax[0].plot(v,t0, 'o', color='w', ms=15, mew=3, lw=10)
    
    if len(Veloo) > 0 and len(Timee) > 0:
        ax[0].plot(Veloo,Timee, '*', color='w', mec='k', mew=2, ms=18, label='Picked location')
        if len(Veloo) > 1 and len(Timee) > 1:
            interpolation = interp1d(Timee, Veloo, kind='linear', fill_value='extrapolate')
            ax[0].plot(interpolation(time), time, color='yellow', label='Interpolated Velocity')
        
    ax[0].vlines(v, np.min(time), np.max(time), colors='w')
    ax[0].hlines(t0, np.min(vel_range), np.max(vel_range), colors='w')
    ax[0].set_xlim(np.min(vel_range), np.max(vel_range))
    ax[0].set_ylim(np.max(time), np.min(time))
    if len(Veloo) == 0 and len(Timee) == 0:
        ax[0].plot([],[],alpha=0,label=' ')
    ax[0].legend(loc='upper right')
    
    plot_seismogram(ax[1], Shift_NMO, offset, time, plot_option=1)
    ax[1].set_title('CMP Gather', fontsize=14, fontweight='bold')
    ax[1].set_ylabel('')

    plot_seismogram(ax[2], Shift_NMO, offset, time, plot_option=3)
    ax[2].set_title('CMP Stack', fontsize=14, fontweight='bold')

    plt.tight_layout()

#this for normalizing np array
def NormalizeData(data):
    data = data.astype(float)
    min_val = np.min(data)
    max_val = np.max(data)
    data = (data - min_val) / (max_val - min_val)
    data = data * 2 - 1
    return data
#-----------------------------------------------------------------------------------------------------------------------------#
# this transfare fruncation
def CrossCoherence (dataRes,dataSor):
  FFTtestRes=fft(np.copy(dataRes))
  FFTtestSoc=fft(np.copy(dataSor))
  upcoherence=(FFTtestRes/np.abs(FFTtestRes))
  lowercoherence=(FFTtestSoc/np.abs(FFTtestSoc))
  coherence=ifft(upcoherence/lowercoherence)
  return coherence
#-----------------------------------------------------------------------------------------------------------------------------#
def CrossCorrelate (dataRes,dataSor):
  testRes=np.copy(dataRes)
  testSoc=np.copy(dataSor)

  corr = correlate( testRes, testSoc,mode='same', method='fft')
  return (corr)
#-----------------------------------------------------------------------------------------------------------------------------#
def CrosscoherenceAuto (dataRes,dataSor):
  testRes=np.copy(dataRes)
  testSoc=np.copy(dataSor)

  corr = coherence(testRes, testSoc, nperseg=7498)
  corr = ifft(corr)
  return (corr)
#-----------------------------------------------------------------------------------------------------------------------------#
# these two fruncations for stacking
def Gup (Gk):
    Gup=(1/np.var(Gk[-1500:-200:]))*Gk
    return Gup
def Wdown (Gj):
    Wdown=(1/np.var(Gj[-1500:-200:]))
    return Wdown
#-----------------------------------------------------------------------------------------------------------------------------#
# these two fruncations for stacking using RMS
def GupRMS (Gk):
    Gup=(np.sqrt(np.mean(Gk[750:1500:]**2))/np.sqrt(np.mean(Gk[-500:-200:]**2))*Gk)
    return Gup
def WdownRMS (Gj):
    Wdown=(np.sqrt(np.mean(Gj[750:1500:]**2))/np.sqrt(np.mean(Gj[-500:-200:]**2)))
    return Wdown
#-----------------------------------------------------------------------------------------------------------------------------#
# this for low band filter cut and band pass filter for nummpy list
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
#-----------------------------------------------------------------------------------------------------------------------------#
# this for coommen noise
def CaluComNoise (data,first,last,channels):
    Sum =np.zeros([7500])
    for x in range (first,last):
        Sum = Sum + np.copy(data[x])

    AvgSum = np.divide(Sum, (last-first))

    for x in range(channels):
        data[x]=np.subtract(data[x],AvgSum)
    return data
#-----------------------------------------------------------------------------------------------------------------------------#
# this for coommen noise for speacl usges
def CaluComNoiseSpe (data,first,last,channels):
    Sum =np.zeros([250])
    for x in range (first,last):
        npdata=np.copy(data[x])
        npdata=np.roll(npdata, 125)
        Sum = Sum + npdata[0:250:]

    AvgSum1 = np.divide(Sum, (last-first))
    AvgSum2 = np.zeros([7250])
    AvgSum = np.append(AvgSum1,AvgSum2)
    AvgSum=np.roll(AvgSum, -125)

    for x in range(channels):
        data[x,]=np.subtract(data[x],AvgSum)
    return data
#-----------------------------------------------------------------------------------------------------------------------------#
# this for calcute time of run
def CalTime():
    now = datetime.now()
    min = now.minute
    sec = now.second
    time = int(min)*60+int(sec)
    return time
#-----------------------------------------------------------------------------------------------------------------------------#
# this for calcute the distance between trances
def distanceOffset (trc1,trc2):
  lat1=trc1.stats.sac['stla']
  lat2=trc2.stats.sac['stla']

  lon1=trc1.stats.sac['stlo']
  lon2=trc2.stats.sac['stlo']

  dist = mpu.haversine_distance((lat1, lon1), (lat2, lon2))*1000
  dist = int(dist)
  return dist

#-----------------------------------------------------------------------------------------------------------------------------#
# this for nomrlize based on speed
def NormSpeed(data,o,stationDistnace,velocty,j,windowNoR):
  velocty = (velocty)/1000*2
  windowNoR = int(windowNoR*1000/2+j)
  data[int(abs(stationDistnace[o])/velocty)+j:int(abs(stationDistnace[o]/velocty))+windowNoR:
     ]=NormalizeData(data[int(abs(stationDistnace[o])/velocty)+j:int(abs(stationDistnace[o]/velocty))+windowNoR:])
  return data
#-----------------------------------------------------------------------------------------------------------------------------#
def readH5(dirct):
  st=Stream()
  f = h5py.File(dirct, 'r')
  dataset = f['Acquisition/Raw[0]/RawData']
  DateTime = f['Acquisition/Raw[0]/RawDataTime']
  DateTimeUTC=(UTCDateTime(np.copy(DateTime[0])/1000000))
  d=np.copy(dataset)
  for x in range(1216):
    a=d[ :,x]
    TestTrace=Trace(data=a)
    TestTrace.stats.sampling_rate=1000
    TestTrace.stats.station=str(x)
    TestTrace.stats.starttime=DateTimeUTC
    st+=TestTrace
  return st
#-----------------------------------------------------------------------------------------------------------------------------#
def readAsccii(dirctor,shfit):
  f = open(dirctor, 'r') # 'r' = read
  ASCdata=[]
  for line in f:
    ASCdata=np.append(ASCdata,float(line)) # note, coma erases the "cartridge return"
  f.close()
  TASCdata = Trace(data=ASCdata)

  Time = dirctor[-15::].replace("-", "T")
  TASCdata.stats.starttime=UTCDateTime(Time)-UTCDateTime(9*60*60)+shfit
  TASCdata.stats.sampling_rate=1000
  return TASCdata
#-----------------------------------------------------------------------------------------------------------------------------#
def apply_minimum_phase(trace):
    # Compute the analytic signal using Hilbert transform
    analytic_signal = hilbert(trace)

    # Take the natural logarithm of the absolute value of the analytic signal
    log_amplitude = np.log(np.abs(analytic_signal))

    # Apply inverse Fourier transform to obtain the minimum phase trace
    minimum_phase_trace = np.real(np.fft.ifft(log_amplitude))

    # Resample the minimum phase trace to the original length
    minimum_phase_trace = resample(minimum_phase_trace, len(trace))

    return minimum_phase_trace
#-----------------------------------------------------------------------------------------------------------------------------#
def calculate_snr(trace, signal_window_start=0, signal_window_end=41):
    """Calculates the signal-to-noise ratio of a trace.

    Args:
        trace: An obspy Trace object.
        signal_window_start: Start time of the signal window in seconds.
        signal_window_end: End time of the signal window in seconds.

    Returns:
        The SNR of the trace.
    """

    signal = trace.slice(starttime=trace.stats.starttime + signal_window_start,
                        endtime=trace.stats.starttime + signal_window_end)
    noise = trace.slice(starttime=trace.stats.endtime - signal_window_end,
                        endtime=trace.stats.endtime)

    signal_variance = np.var(signal.data)
    noise_variance = np.var(noise.data)

    if noise_variance == 0:
        return float('inf')  # Handle cases where noise variance is zero
    else:
        return signal_variance / noise_variance
#-----------------------------------------------------------------------------------------------------------------------------#
def extract_seismic_data(file_path):
    with open(file_path, 'r', encoding='shift_jis') as file:
        lines = file.readlines()

    # Find the start of the data section
    data_section_start = None
    for i, line in enumerate(lines):
        if '[DATA]' in line:
            data_section_start = i + 1
            break

    # If data section exists, extract the data
    if data_section_start:
        data_lines = lines[data_section_start:]
        data = []
        for line in data_lines:
            values = line.split(',')
            if len(values) >= 4:
                try:
                    # Extract the first column (or other relevant data)
                    data.append(float(values[0]))  # Modify index if needed
                except ValueError:
                    continue  # Skip invalid lines
        return np.array(data)

    return None
    #-----------------------------------------------------------------------------------------------------------------------------#
def extract_seismic_Time(file_path):
    with open(file_path, 'r', encoding='shift_jis') as file:
        lines = file.readlines()
        DateTime= lines[1][9:]
        
    return (UTCDateTime(DateTime))

#-----------------------------------------------------------------------------------------------------------------------------#

def weighted_stack(stream):
    """Stacks traces in a stream, weighting them by their SNR
       and grouping by trace ID.
    """

    stacked_traces = {}
    for tr in stream:
        trace_id = tr.id  # Get the trace ID (e.g., 'station.location')
        if trace_id not in stacked_traces:
            stacked_traces[trace_id] = []
        stacked_traces[trace_id].append(tr)

    result_stream = Stream()
    for trace_id, traces in stacked_traces.items():
        # Calculate weights for traces with the same ID
        weights = np.array([calculate_snr(tr) for tr in traces])
        normalized_weights = weights / np.sum(weights)

        # Perform weighted stacking for this group
        stacked_trace = np.zeros(len(traces[0].data))
        for i, tr in enumerate(traces):
            stacked_trace += tr.data * normalized_weights[i]

        # Create a new Trace object for the stacked data
        stacked_stream = traces[0].copy()
        stacked_stream.data = stacked_trace
        result_stream += stacked_stream

    return result_stream
#-----------------------------------------------------------------------------------------------------------------------------#
def apply_agc(stream, window_length=0.5):
    """
    Apply Automatic Gain Control (AGC) to each trace in the ObsPy stream.
    
    Parameters:
    - stream: obspy.Stream, the seismic data stream to process.
    - window_length: float, the AGC window length in seconds.
    
    Returns:
    - stream: obspy.Stream, the AGC-processed stream.
    """
    for trace in stream:
        data = trace.data.astype(np.float64)  # Ensure data is in float format
        sample_rate = trace.stats.sampling_rate
        samples_per_window = int(window_length * sample_rate)
        
        # Calculate the envelope manually using the Hilbert transform
        data_envelope = np.abs(hilbert(data))
        
        # Compute a moving average over the envelope
        smoothed_envelope = np.convolve(data_envelope, np.ones(samples_per_window), mode='same') / samples_per_window
        
        # Apply AGC by dividing by the smoothed envelope
        data_agc = data / (smoothed_envelope + 1e-12)  # Add a small constant to avoid division by zero
        trace.data = data_agc
    
    return stream


