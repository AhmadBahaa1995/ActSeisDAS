"""
EdgeDAS: An Automated Onsite Processing Framework for Distributed Acoustic Sensing (DAS) Data

This script provides a comprehensive framework for the automated processing of Distributed
Acoustic Sensing (DAS) data. It is designed to run on an edge computing device located
onsite, enabling near-real-time analysis of seismic data.

The framework includes the following key features:
- A graphical user interface (GUI) built with Tkinter for easy configuration.
- Automated file discovery and processing of raw binary DAS data.
- Parallel processing capabilities using multiprocessing to leverage multi-core CPUs
  for efficient cross-correlation.
- Cross-correlation of continuous DAS data against a reference source signal to
  generate virtual shot gathers.
- Options for both standard and binary (one-bit) cross-correlation.
- Stacking of shot gathers over user-defined time periods (e.g., hours) to improve
  the signal-to-noise ratio.
- Automated generation and saving of seismogram plots for quality control.
- A scheduling feature to run the processing pipeline automatically at a specified
  time each day.
- Persistent settings, allowing the application to remember the last used
  configuration.
- **Robust, Process-Safe Logging**: Uses Python's standard logging module to prevent
  file access errors during parallel processing.
- **Performance Monitoring**: Detailed logging of time, CPU, and memory usage for
  each processing stage into a structured CSV file for analysis.

Workflow:
1. The user configures processing parameters via the GUI (e.g., data folder,
   shot interval, stacking duration).
2. Upon clicking "Run," the application starts an initial processing run and
   schedules subsequent daily runs in a background thread to keep the GUI responsive.
3. The script identifies unprocessed days by comparing the data folders against a log file.
4. For each unprocessed day, it distributes the workload across multiple CPU cores. Each
   core processes one hour of data, recording performance metrics for each step.
5. Worker processes write their output files directly to the final `_Processed_Data`
   directory.
6. The script aggregates performance data and writes it to `performance_log.csv`.
7. The script then stacks the hourly shot gathers from the final destination.
8. Finally, it generates a seismogram plot, saves the final data, logs the
   performance of these final steps, and marks the day as complete.
"""

import logging
import os
import struct
import numpy as np
import schedule
import time
import matplotlib.pyplot as plt
import gc
import sys
import psutil
import tkinter as tk
import threading
import re
import csv
import tempfile
import shutil
import logging.handlers

from tkinter import filedialog, messagebox, OptionMenu, StringVar, BooleanVar, Checkbutton, Radiobutton, Entry, DISABLED, NORMAL
from scipy.signal import butter, lfilter, correlate, coherence, minimum_phase, hilbert, find_peaks
from scipy.fftpack import fft, fftfreq
from obspy import UTCDateTime, read, Stream, Trace
from obspy.signal import cross_correlation
from multiprocessing import Pool, cpu_count, Manager, Queue
from tqdm import tqdm
from datetime import datetime
from functools import partial

# --- Global Setup ---
# Get the directory where the script is located
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.getcwd() # Fallback for interactive environments

# Define a LOCAL log directory to avoid network issues and for cleaner organization
local_log_dir = os.path.join(script_dir, '_log')
os.makedirs(local_log_dir, exist_ok=True)
log_file_path = os.path.join(local_log_dir, 'processing_log.txt')

# --- Multiprocessing-safe Logging Setup ---
def setup_logging_listener(log_queue):
    """
    Configures and starts a listener in the main process.
    The listener pulls log records from a queue and sends them to the configured handlers.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(processName)s - %(levelname)s - %(message)s')

    file_handler = logging.FileHandler(log_file_path, mode='a', encoding='utf-8')
    file_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)

    queue_listener = logging.handlers.QueueListener(log_queue, file_handler, console_handler)
    queue_listener.start()
    return queue_listener

def worker_logging_init(log_queue):
    """
    Configures logging for each worker process.
    It removes any existing handlers and adds a QueueHandler to send all logs
    back to the main process's listener queue.
    """
    queue_handler = logging.handlers.QueueHandler(log_queue)
    root_logger = logging.getLogger()
    root_logger.handlers.clear()
    root_logger.addHandler(queue_handler)
    root_logger.setLevel(logging.INFO)


# --------------------------------------------------------------------------------------------------------- #
def find_list_difference(list1, list2):
    """
    Finds elements that are in list1 but not in list2.

    Args:
        list1 (list): The primary list.
        list2 (list): The list of elements to exclude.

    Returns:
        list: A list containing elements unique to list1.
    """
    return list(set(list1) - set(list2))

# --------------------------------------------------------------------------------------------------------- #
def plot_seismogram(ax, data, offset, time, plot_option=1, cmap='seismic', title=None):
    """
    A versatile function to plot seismic data in various formats.
    """
    # --- Nested plotting functions for organization ---
    def wiggle(ax, data, time):
        """Plots data as wiggle traces."""
        for i in range(data.shape[1]):
            trace_norm = data[:, i] / (np.max(np.abs(data[:, i])) + 1e-9)
            ax.plot(trace_norm + i, time, linewidth=0.3, c='k')
            ax.fill_betweenx(time, trace_norm + i, i, where=(trace_norm + i >= i), color='k', lw=0)
        ax.invert_yaxis()
        ax.set_ylabel('Time [s]')
        ax.set_xlabel('Trace Number')

    def color(ax, data, time, offset):
        """Plots data as a color image."""
        clim = np.percentile(np.abs(data), 98)
        ax.imshow(data, cmap=cmap,
                  extent=[np.min(offset), np.max(offset), np.max(time), np.min(time)],
                  aspect='auto', vmin=-clim, vmax=clim)
        ax.set_ylabel('Time [s]')
        ax.set_xlabel('Trace Number')

    def stack_seismogram(ax, data, time):
        """Plots the stacked (summed) trace of all channels."""
        stacked_trace = data.sum(axis=1)
        ax.plot(stacked_trace, time, lw=0.5, color='k')
        ax.fill_betweenx(time, stacked_trace, 0, where=(stacked_trace > 0), color='k', lw=0)
        ax.set_xlabel('Stacked Amplitude')
        ax.set_ylabel('Time [s]')
        ax.invert_yaxis()

    # --- Main plotting logic ---
    plt.rcParams["figure.figsize"] = [12, 6]
    if plot_option == 0 or plot_option == 'wiggle':
        wiggle(ax, data, time)
    elif plot_option == 1 or plot_option == 'color':
        color(ax, data, time, offset)
    elif plot_option == 2 or plot_option == 'both':
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(15, 7))
        wiggle(ax1, data, time)
        color(ax2, data, time, offset)
        if title:
            fig.suptitle(title)
        return
    elif plot_option == 3 or plot_option == 'stacked':
        stack_seismogram(ax, data, time)

    if title:
        ax.set_title(title)

# --------------------------------------------------------------------------------------------------------- #
def log_performance_metric(log_path, metric):
    """
    Appends a single performance metric dictionary to a CSV file.
    """
    header = [
        "timestamp", "day", "hour", "stage", "duration_seconds",
        "cpu_percent", "mem_rss_mb_start", "mem_rss_mb_end", "mem_delta_mb"
    ]
    file_exists = os.path.isfile(log_path)
    try:
        with open(log_path, 'a', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=header)
            if not file_exists:
                writer.writeheader()
            writer.writerow(metric)
    except Exception as e:
        logging.error(f"Error writing to performance log: {e}")

# --------------------------------------------------------------------------------------------------------- #
def plot_CC(merged_data, target_day, StartTime, EndTime, output_base_path):
    """
    Plots and saves the final stacked cross-correlation seismogram.
    """
    try:
        t_start = time.perf_counter()
        process = psutil.Process(os.getpid())
        mem_start = process.memory_info().rss / (1024 * 1024)

        start_str = StartTime.strftime('%Y-%m-%dT%H-%M-%S')
        end_str = EndTime.strftime('%H-%M-%S')
        title = f"Stacked CC: {target_day} from {start_str} to {end_str}"
        safe_title = title.replace(":", "-").replace(" ", "_")

        logging.info(f"Creating plot: {title}")
        merged_data.normalize()
        fig, ax = plt.subplots(figsize=(12, 6))
        plot_data = np.array([tr.data for tr in merged_data]).T
        sr = merged_data[0].stats.sampling_rate
        npts = merged_data[0].stats.npts
        time_axis = np.arange(-0.5, npts / sr - 0.5, 1 / sr)
        offset_axis = np.arange(len(merged_data))

        plot_seismogram(ax, plot_data, offset=offset_axis, time=time_axis, cmap="seismic", plot_option=1, title=title)
        
        output_directory = os.path.join(output_base_path, "Seismogram_Plots")
        os.makedirs(output_directory, exist_ok=True)
        save_path = os.path.join(output_directory, safe_title + ".png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        logging.info(f"Saved seismogram plot at {save_path}")

        cc_output_directory = os.path.join(output_base_path, "Stacked_Cross_Correlation")
        os.makedirs(cc_output_directory, exist_ok=True)
        file_path = os.path.join(cc_output_directory, f"{safe_title}.mseed")
        merged_data.write(file_path, format="MSEED")
        logging.info(f"Saved stacked CC data at {file_path}")
        
        duration = time.perf_counter() - t_start
        mem_end = process.memory_info().rss / (1024 * 1024)
        perf_log_path = os.path.join(output_base_path, 'performance_log.csv')
        metric = {
            "timestamp": datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "day": target_day,
            "hour": "N/A", "stage": "Plotting_and_Saving_Stacked", "duration_seconds": round(duration, 4),
            "cpu_percent": psutil.cpu_percent(), "mem_rss_mb_start": round(mem_start, 2),
            "mem_rss_mb_end": round(mem_end, 2), "mem_delta_mb": round(mem_end - mem_start, 2)
        }
        log_performance_metric(perf_log_path, metric)

    except Exception as e:
        logging.error(f"Error creating plot or saving data for {title}: {e}")

    finally:
        log_path = os.path.join(local_log_dir, "donedays.txt")
        with open(log_path, "a", encoding='utf-8') as log_file:
            log_file.write(f"{target_day}, finished at: {datetime.now()}\n")
        logging.info(f"Finished processing for day: {target_day}")

# --------------------------------------------------------------------------------------------------------- #
def stacking(stacked_hours, available_Hr_Path, target_day, output_base_path):
    """
    Reads hourly cross-correlation files and stacks them over a specified period.
    """
    t_start = time.perf_counter()
    process = psutil.Process(os.getpid())
    mem_start = process.memory_info().rss / (1024 * 1024)

    available_hr_folders = sorted(os.listdir(available_Hr_Path))
    data_container = Stream()
    read_hr_counter = 0
    shots_generated = 0

    logging.info(f"Starting stacking process for data in {available_Hr_Path}")
    for hr_folder in tqdm(available_hr_folders, desc="Stacking hours"):
        hr_path = os.path.join(available_Hr_Path, hr_folder)
        if not os.path.isdir(hr_path):
            continue
        for data_file in sorted(os.listdir(hr_path)):
            file_path = os.path.join(hr_path, data_file)
            if os.path.isfile(file_path) and data_file.endswith(".mseed"):
                try:
                    data_container += read(file_path)
                except Exception as e:
                    logging.warning(f"Could not read file {file_path}. Error: {e}")
        read_hr_counter += 1

        if len(data_container) > 0 and read_hr_counter % stacked_hours == 0:
            logging.info(f"Stacking data for a {read_hr_counter}-hour block...")
            StartTime = data_container[0].stats.starttime
            EndTime = data_container[-1].stats.endtime
            
            data_container.stack('id')
            shots_generated += 1
            plot_CC(data_container, target_day, StartTime, EndTime, output_base_path)
            del data_container; gc.collect()
            data_container = Stream()
            read_hr_counter = 0
    
    if len(data_container) > 0:
        logging.info(f"Stacking remaining {read_hr_counter} hour(s) of data...")
        StartTime = data_container[0].stats.starttime
        EndTime = data_container[-1].stats.endtime
        data_container.stack('id')
        shots_generated += 1
        plot_CC(data_container, target_day, StartTime, EndTime, output_base_path)
        del data_container; gc.collect()
        
    duration = time.perf_counter() - t_start
    mem_end = process.memory_info().rss / (1024 * 1024)
    perf_log_path = os.path.join(output_base_path, 'performance_log.csv')
    metric = {
        "timestamp": datetime.now().strftime('%Y-%m-%d %H:%M:%S'), "day": target_day,
        "hour": "N/A", "stage": "Stacking_All_Hours", "duration_seconds": round(duration, 4),
        "cpu_percent": psutil.cpu_percent(), "mem_rss_mb_start": round(mem_start, 2),
        "mem_rss_mb_end": round(mem_end, 2), "mem_delta_mb": round(mem_end - mem_start, 2)
    }
    log_performance_metric(perf_log_path, metric)
    
    logging.info(f"Stacking complete for day: {target_day}. Generated {shots_generated} stacked seismic shot gathers.")
    return None

# --------------------------------------------------------------------------------------------------------- #
def process_channel_chunk(args):
    """
    Worker function to process one hour of raw binary data using a sliding window approach.
    It writes output directly to the final destination folder.
    """
    hour_folder, shotting_interval, day_folder, shift_seconds, source_pass_or_mode, use_binary_cc, output_base_path, source_channel, min_freq, max_freq = args
    
    day_str = os.path.basename(day_folder)
    hour_str = os.path.basename(hour_folder)
    save_path = os.path.join(output_base_path, 'Cross_Correlation_folder', day_str, hour_str)
    source_save_dir = os.path.join(output_base_path, 'Source_Signals', day_str, hour_str)
    os.makedirs(save_path, exist_ok=True)
    os.makedirs(source_save_dir, exist_ok=True)
    
    metrics = []

    try:
        file_list = sorted(os.listdir(hour_folder))

        # --- File Skipping Logic ---
        process_list = file_list[1:-1] 
        files_to_process = []
        if shotting_interval > 60:
            step = max(1, int(round(shotting_interval / 60.0)))
    
            files_to_process = process_list[::step]
        else:
            files_to_process = process_list

        for file_name in tqdm(files_to_process, desc=f'Processing {hour_str}', leave=False):
            file_path = os.path.join(hour_folder, file_name)
            if not os.path.isfile(file_path): continue

            try:
                # --- Stage: I/O Read and Data Preparation ---
                clean_file_name = file_name.replace("_EXTTRIG", "")
                dt = datetime.strptime(clean_file_name, "%Y%m%d%H%M%S%f")
                with open(file_path, 'br') as f:
                    data = []
                    num_channels_expected = 1200
                    num_samples_per_channel = 0

                    for _ in range(num_channels_expected):
                        header = f.read(8)
                        if not header: break
                        num_sections, num_samples = struct.unpack('<ii', header)
                        if num_samples_per_channel == 0: num_samples_per_channel = num_samples
                        
                        total_samples = num_sections * num_samples
                        channel_raw_data = f.read(total_samples * 4)
                        data.extend(list(struct.unpack(f'<{total_samples}f', channel_raw_data)))

                if not data: continue

                # --- Source Preparation ---
                source_pass = None
                if source_pass_or_mode == 'extract':
                    raw_source_data = data[source_channel::num_samples_per_channel]
                    source_trace = Trace(data=np.array(raw_source_data)[::2])
                    source_trace.stats.sampling_rate = 500
                    source_trace.stats.starttime = UTCDateTime(dt)
                    source_trace.detrend()
                    source_trace.filter('bandpass', freqmin=min_freq, freqmax=max_freq, corners=8)
                    source_trace.stats.station = f"{source_channel:04d}"
                    source_pass = source_trace
                else:
                    source_pass = source_pass_or_mode[0]

                # --- Main Processing Loop ---
                all_channels_cc_stream = Stream()
                channels_data = [data[i::num_samples_per_channel] for i in range(num_samples_per_channel)]

                
                total_windows_processed = 0

                for i, channel_data in enumerate(channels_data):
                    trace = Trace(data=np.array(channel_data)[::2])
                    trace.stats.station = f"{i:04d}"
                    trace.stats.sampling_rate = 500
                    trace.stats.starttime = UTCDateTime(dt)
                    trace.detrend('linear')
                    trace.filter('bandpass', freqmin=min_freq, freqmax=max_freq, corners=8)


                    for window in trace.slide(window_length=55.996, step=shotting_interval,
                                              offset=shift_seconds,
                                              include_partial_windows=False, nearest_sample=True):
                        total_windows_processed += 1
                        current_source = source_pass.copy().trim(starttime=window.stats.starttime, endtime=window.stats.endtime)
                        
                        if len(current_source.data) == len(window.data):
                            if use_binary_cc:
                                window_rms = np.sqrt(np.mean(window.data**2))
                                window_thresholded = np.where(np.abs(window.data) > window_rms, np.sign(window.data).astype(int), 0)
                                source_pass_rms = np.sqrt(np.mean(current_source.data**2))
                                source_pass_thresholded = np.where(np.abs(current_source.data) > source_pass_rms, np.sign(current_source.data).astype(int), 0)
                                cc_data = cross_correlation.correlate(window_thresholded, source_pass_thresholded, 2500)
                            else:
                                cc_data = cross_correlation.correlate(window, current_source, 2500)
                            
                            if cc_data is not None and cc_data.size > 0:
                                cc_trace = Trace(data=cc_data[2250:])
                                cc_trace.stats.station = window.stats.station
                                cc_trace.stats.starttime = window.stats.starttime
                                cc_trace.stats.sampling_rate = window.stats.sampling_rate
                                all_channels_cc_stream += cc_trace
                
                if len(all_channels_cc_stream) > 0:
                    formatted_starttime = str(all_channels_cc_stream[0].stats.starttime).replace(':', '-').split('.')[0]
                    output_filename = f'{formatted_starttime}.mseed'
                    all_channels_cc_stream.write(os.path.join(save_path, output_filename), format='MSEED')
                
                    if source_pass_or_mode == 'extract':
                        source_filename = f"{clean_file_name}_source_ch{source_channel}.mseed"
                        source_pass.write(os.path.join(source_save_dir, source_filename), format="MSEED")

                del data, all_channels_cc_stream, channels_data
                gc.collect()

            except Exception as e:
                logging.error(f"Error processing file {file_name}: {e}")

    except Exception as e:
        logging.error(f"Error processing hour folder {hour_folder}: {e}")

    return metrics

# --------------------------------------------------------------------------------------------------------- #
def distribute_work(day_folder, shift_seconds, shotting_interval, target_day, stacked_hours, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq, log_queue):
    """
    Distributes the processing of one day's data across multiple CPU cores.
    """
    source_pass_or_mode = 'extract'
    if source_mode == 'file':
        source_file = os.path.join(local_log_dir, 'source_pass.mseed')
        if not os.path.exists(source_file):
            messagebox.showerror("Error", f"Source file not found:\n{source_file}")
            return
        source_pass = read(source_file)
        source_pass_or_mode = source_pass

    hour_list = sorted([h for h in os.listdir(day_folder) if os.path.isdir(os.path.join(day_folder, h)) and h.isdigit() and len(h) == 2])
    if not hour_list:
        logging.warning(f"No valid hour folders found in {day_folder}. Skipping.")
        return
    
    cores = min(cores, len(hour_list))
    logging.info(f"Using {cores} cores for parallel processing of {len(hour_list)} hours.")

    final_output_base_path = os.path.join(script_dir, '_Processed_Data')
    os.makedirs(final_output_base_path, exist_ok=True)
    perf_log_path = os.path.join(final_output_base_path, 'performance_log.csv')

    args = [(os.path.join(day_folder, hour), shotting_interval, day_folder, shift_seconds, source_pass_or_mode, use_binary_cc, final_output_base_path, source_channel, min_freq, max_freq) for hour in hour_list]
    
    with Pool(processes=cores, initializer=worker_logging_init, initargs=(log_queue,)) as pool:
        all_metrics = pool.map(process_channel_chunk, args)
    
    logging.info("Parallel processing complete.")

    logging.info("Aggregating performance metrics from all workers...")
    flat_metrics = [item for sublist in all_metrics for item in sublist]
    for metric in flat_metrics:
        log_performance_metric(perf_log_path, metric)
    logging.info(f"Performance metrics saved to {perf_log_path}")

    available_Hr_Path = os.path.join(final_output_base_path, 'Cross_Correlation_folder', os.path.basename(day_folder))
    if os.path.exists(available_Hr_Path):
        available_hr = sorted(os.listdir(available_Hr_Path))
        hours_to_stack = stacked_hours
        if len(available_hr) < stacked_hours:
            logging.warning(f"Only {len(available_hr)} hours of data available for stacking.")
            hours_to_stack = len(available_hr)
        
        if hours_to_stack > 0:
            stacking(hours_to_stack, available_Hr_Path, target_day, final_output_base_path)
        else:
            logging.info("No data available to stack.")
    else:
        logging.warning(f"No cross-correlation data found for {target_day}. Skipping stacking.")

    return None

# --------------------------------------------------------------------------------------------------------- #
def process_seismic_data(shift_seconds, shotting_interval, target_recorded, stacked_hours, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq, log_queue):
    """
    Main function to orchestrate the entire data processing workflow.
    """
    logging.info("--- Starting Seismic Data Processing Workflow ---")
    
    day_list = sorted([d for d in os.listdir(target_recorded) if os.path.isdir(os.path.join(target_recorded, d)) and d.isdigit() and len(d) == 8])
    log_path = os.path.join(local_log_dir, 'donedays.txt')
    if not os.path.exists(log_path): open(log_path, 'w', encoding='utf-8').close()

    with open(log_path, 'r', encoding='utf-8') as f:
        done_days = [line.strip().split(',')[0] for line in f if line.strip()]

    target_days = sorted(find_list_difference(day_list, done_days))
    if not target_days:
        logging.info("No new days to process. All available data has been processed.")
        return

    logging.info(f"Found {len(target_days)} new day(s) to process:  from {target_days[0]}, to {target_days[-1]}")

    for target_day in tqdm(target_days, desc="Processing Days"):
        logging.info(f"--- Starting processing for day: {target_day} ---")
        day_folder = os.path.join(target_recorded, target_day)
        distribute_work(day_folder, shift_seconds, shotting_interval, target_day, stacked_hours, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq, log_queue)

    logging.info("--- All new seismic data has been processed. ---")


# --------------------------------------------------------------------------------------------------------- #
def schedule_daily_processing(run_time, shift_seconds, shotting_interval, target_folder, hours_to_stack, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq, log_queue):
    """
    Schedules the main processing function to run daily at a specified time.
    """
    schedule.every().day.at(run_time).do(process_seismic_data, shift_seconds, shotting_interval, target_folder, hours_to_stack, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq, log_queue)
    logging.info(f"Processing has been scheduled to run daily at {run_time}.")
    while True:
        schedule.run_pending()
        time.sleep(1)

# --------------------------------------------------------------------------------------------------------- #
def save_settings(shift_seconds, shotting_interval, target_folder, hours_to_stack, run_time, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq):
    """Saves the current GUI settings to a file for persistence."""
    try:
        with open(SETTINGS_FILE, 'w', encoding='utf-8') as f:
            f.write(f"shift_seconds={shift_seconds}\n")
            f.write(f"shotting_interval={shotting_interval}\n")
            f.write(f"target_folder={target_folder}\n")
            f.write(f"hours_to_stack={hours_to_stack}\n")
            f.write(f"run_time={run_time}\n")
            f.write(f"use_binary_cc={use_binary_cc}\n")
            f.write(f"cores={cores}\n")
            f.write(f"source_mode={source_mode}\n")
            f.write(f"source_channel={source_channel}\n")
            f.write(f"min_freq={min_freq}\n")
            f.write(f"max_freq={max_freq}\n")
        logging.info("Processing settings saved.")
    except Exception as e:
        logging.error(f"Error saving settings: {e}")

# --------------------------------------------------------------------------------------------------------- #
def load_settings():
    """Loads processing settings from the file to populate the GUI."""
    settings = {}
    if os.path.exists(SETTINGS_FILE):
        try:
            with open(SETTINGS_FILE, 'r', encoding='utf-8') as f:
                for line in f:
                    if "=" in line:
                        key, value = line.strip().split("=", 1)
                        settings[key] = value
            logging.info("Processing settings loaded from file.")
        except Exception as e:
            logging.error(f"Error loading settings: {e}")
    return settings

# --------------------------------------------------------------------------------------------------------- #
# --- GUI Application ---
# --------------------------------------------------------------------------------------------------------- #
if __name__ == "__main__":
    SETTINGS_FILE = os.path.join(local_log_dir, "processing_settings.txt")
    
    # Create a manager queue for logging from multiple processes
    log_queue = Manager().Queue()
    # Start the listener that will handle logs from all processes
    queue_listener = setup_logging_listener(log_queue)


    def start_processing_thread():
        """
        Wrapper function to get GUI values and start the processing and scheduling
        in a background thread to prevent the GUI from freezing.
        """
        try:
            # --- Get all values from GUI ---
            shift_seconds = int(shift_seconds_entry.get())
            shotting_interval = int(shotting_interval_entry.get())
            target_folder = target_folder_entry.get()
            hours_to_stack = int(hours_stack_var.get() or 24)
            selected_hour = hour_var.get()
            selected_minute = minute_var.get()
            run_time = f"{selected_hour}:{selected_minute}"
            cores = int(cores_entry.get() or 24)
            use_binary_cc = binary_cc_var.get()
            source_mode = source_mode_var.get()
            source_channel = int(source_channel_entry.get()) if source_mode == 'extract' else 0
            min_freq = float(min_freq_entry.get())
            max_freq = float(max_freq_entry.get())


            if not os.path.isdir(target_folder):
                messagebox.showerror("Error", "Target folder path is invalid.")
                return

            save_settings(shift_seconds, shotting_interval, target_folder, hours_to_stack, run_time, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq)
            messagebox.showinfo("Processing Started", "Seismic data processing has started in the background.\nThe application will also run automatically every day.")

            def background_task():
                process_seismic_data(shift_seconds, shotting_interval, target_folder, hours_to_stack, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq, log_queue)
                schedule_daily_processing(run_time, shift_seconds, shotting_interval, target_folder, hours_to_stack, use_binary_cc, cores, source_mode, source_channel, min_freq, max_freq, log_queue)

            processing_thread = threading.Thread(target=background_task, daemon=True)
            processing_thread.start()

        except ValueError:
            messagebox.showerror("Error", "Please enter valid numbers for all numeric fields (Shift, Interval, Cores, Frequencies, etc.).")
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

    def browse_folder():
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            target_folder_entry.delete(0, tk.END)
            target_folder_entry.insert(0, folder_selected)

    def toggle_source_channel_entry():
        source_channel_entry.config(state=NORMAL if source_mode_var.get() == "extract" else DISABLED)

    # --- Build the GUI Window ---
    window = tk.Tk()
    window.title("EdgeDAS: Onsite Processing Framework")

    # --- Configuration Frame ---
    config_frame = tk.LabelFrame(window, text="Processing Configuration", padx=10, pady=10)
    config_frame.pack(padx=10, pady=10, fill="x")

    tk.Label(config_frame, text="Target Data Folder:").grid(row=0, column=0, sticky="w", pady=2)
    target_folder_entry = tk.Entry(config_frame, width=30)
    target_folder_entry.grid(row=0, column=1, padx=5, pady=2, columnspan=3)
    browse_button = tk.Button(config_frame, text="Browse...", command=browse_folder)
    browse_button.grid(row=0, column=4, padx=5, pady=2)

    tk.Label(config_frame, text="Shift Seconds:").grid(row=1, column=0, sticky="w", pady=2)
    shift_seconds_entry = tk.Entry(config_frame, width=10)
    shift_seconds_entry.grid(row=1, column=1, padx=5, pady=2, sticky="w")

    tk.Label(config_frame, text="Shooting Interval (s):").grid(row=2, column=0, sticky="w", pady=2)
    shotting_interval_entry = tk.Entry(config_frame, width=10)
    shotting_interval_entry.grid(row=2, column=1, padx=5, pady=2, sticky="w")

    tk.Label(config_frame, text="Hours to Stack:").grid(row=3, column=0, sticky="w", pady=2)
    hours_stack_var = StringVar(window)
    hours_options = ["1", "2", "4", "6", "8", "12", "24"]
    hours_stack_menu = OptionMenu(config_frame, hours_stack_var, *hours_options)
    hours_stack_menu.grid(row=3, column=1, padx=5, pady=2, sticky="ew")

    # --- Source Function Frame ---
    source_frame = tk.LabelFrame(window, text="Source Function", padx=10, pady=10)
    source_frame.pack(padx=10, pady=10, fill="x")
    
    source_mode_var = StringVar(value="file")
    Radiobutton(source_frame, text="Use Pre-prepared File (source_pass.mseed)", variable=source_mode_var, value="file", command=toggle_source_channel_entry).grid(row=0, column=0, sticky="w", columnspan=2)
    Radiobutton(source_frame, text="Extract from DAS Data (Channel #):", variable=source_mode_var, value="extract", command=toggle_source_channel_entry).grid(row=1, column=0, sticky="w")
    source_channel_entry = Entry(source_frame, width=10)
    source_channel_entry.grid(row=1, column=1, sticky="w", padx=5)

    # --- Advanced Settings Frame ---
    advanced_frame = tk.LabelFrame(window, text="Advanced Settings", padx=10, pady=10)
    advanced_frame.pack(padx=10, pady=10, fill="x")

    tk.Label(advanced_frame, text="Daily Run Time (HH:MM):").grid(row=0, column=0, sticky="w", pady=2)
    hour_var, minute_var = StringVar(window), StringVar(window)
    hours, minutes = [f"{i:02d}" for i in range(24)], [f"{i:02d}" for i in range(60)]
    time_frame = tk.Frame(advanced_frame)
    time_frame.grid(row=0, column=1, pady=2, sticky="w", columnspan=3)
    OptionMenu(time_frame, hour_var, *hours).pack(side="left")
    tk.Label(time_frame, text=":").pack(side="left", padx=2)
    OptionMenu(time_frame, minute_var, *minutes).pack(side="left")

    tk.Label(advanced_frame, text="Cores to Use (24 if blank):").grid(row=1, column=0, sticky="w", pady=2)
    cores_entry = tk.Entry(advanced_frame, width=10)
    cores_entry.grid(row=1, column=1, pady=2, sticky="w")
    
    binary_cc_var = BooleanVar()
    Checkbutton(advanced_frame, text="Use Binary Cross-Correlation", variable=binary_cc_var).grid(row=3, column=0, columnspan=2, sticky="w", pady=5)
    
    # New Filter Frequency entries
    tk.Label(advanced_frame, text="Bandpass Filter:").grid(row=2, column=0, sticky="w", pady=2)
    filter_frame = tk.Frame(advanced_frame)
    filter_frame.grid(row=2, column=1, pady=2, sticky="w", columnspan=3)
    min_freq_entry = tk.Entry(filter_frame, width=6)
    min_freq_entry.pack(side="left")
    tk.Label(filter_frame, text=" - ").pack(side="left", padx=2)
    max_freq_entry = tk.Entry(filter_frame, width=6)
    max_freq_entry.pack(side="left")
    tk.Label(filter_frame, text=" Hz").pack(side="left", padx=2)


    # --- Run Button ---
    run_button = tk.Button(window, text="Run Processing & Schedule Daily", command=start_processing_thread)
    run_button.pack(pady=10, padx=10, ipady=5, fill="x")

    # --- Load settings on startup to populate fields ---
    saved_settings = load_settings()
    if saved_settings:
        shift_seconds_entry.insert(0, saved_settings.get('shift_seconds', '0'))
        shotting_interval_entry.insert(0, saved_settings.get('shotting_interval', '30'))
        target_folder_entry.insert(0, saved_settings.get('target_folder', ''))
        hours_stack_var.set(saved_settings.get('hours_to_stack', '24'))
        run_time = saved_settings.get('run_time', '01:00').split(':')
        hour_var.set(run_time[0])
        minute_var.set(run_time[1])
        cores_entry.insert(0, saved_settings.get('cores', '24'))
        binary_cc_var.set(saved_settings.get('use_binary_cc', 'false').lower() == 'true')
        source_mode_var.set(saved_settings.get('source_mode', 'file'))
        source_channel_entry.insert(0, saved_settings.get('source_channel', '0'))
        min_freq_entry.insert(0, saved_settings.get('min_freq', '15.0'))
        max_freq_entry.insert(0, saved_settings.get('max_freq', '60.0'))
    else: # Set some sensible defaults if no settings file
        shotting_interval_entry.insert(0, '30')
        hours_stack_var.set('24')
        hour_var.set('01')
        minute_var.set('00')
        cores_entry.insert(0, '24')
        source_channel_entry.insert(0, '0')
        min_freq_entry.insert(0, '15.0')
        max_freq_entry.insert(0, '60.0')

    toggle_source_channel_entry() # Set initial state of the entry box
    window.mainloop() # Start the GUI event loop
    queue_listener.stop()


