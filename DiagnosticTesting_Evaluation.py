import io
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, lfilter, freqz
import os.path
import os
import argparse

def encoding_read(filename):
    try:
        df = importdata(filename, 'utf-16 LE')
    except:
        df = importdata(filename, 'utf-8')
    return df

def importdata(filename, encoding):
    # encoding of Pharos Software changed at some point, need to check if code works
    # file = io.open(filename, 'r', encoding='utf-16-le')
    file = io.open(filename, 'r', encoding=encoding)
    list = []
    Cy5 = []
    Fam = []
    Time = []
    cutoff = 0
    print(filename)

    for n, lines in enumerate(file):
        list.append(lines.split())
        if "Melting curve time interval[sec.]" in lines:
            cutoff = n + 40

    # code to extract the channel values into list elements with float values in them
    for n, line in enumerate(list[cutoff:]):
        # make sure that list elements really have TRUE values in them, poor coding
        if len(line) < 10:
            line = list[n - 1]
            print(line)


        Fam.append(float(line[4].replace(',', '.')))
        # Cy5.append(float(line[12].replace(',','.')))
        Time.append(float(line[0].replace(',', '.')))

    df = pd.DataFrame({'TIME': Time, 'FAM': Fam})  # 'NTC-CY5': Cy5})
    # remove baseline, just like in Pharos Software
    df['FAM'] = df['FAM'] - df['FAM'][0]

    return df

def smooth_filt(df, filterlength):
    fl = filterlength
    smooth_signal = np.convolve(df['FAM'], np.ones(fl) / fl, mode='valid')
    # generate new time vector so that x-data and y-data have same length
    time_vec = df.iloc[(int(fl / 2)):-int(fl / 2), 0]
    if len(time_vec) == len(smooth_signal):
        smooth_df = pd.DataFrame({'TIME': time_vec, 'FAM': smooth_signal})
    else:
        time_vec = df.iloc[(int(fl / 2)-1):-int(fl / 2), 0]
        smooth_df = pd.DataFrame({'TIME': time_vec, 'FAM': smooth_signal})
    return smooth_df

def gradient1(df):
    firstder = np.diff(df.iloc[:, 1]) / np.diff(df.iloc[:, 0])
    time_vec = df.iloc[:-1, 0]
    grad_df = pd.DataFrame({'TIME': time_vec, 'FAM': firstder})
    return grad_df

def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def lowpass(data, fs1, cutoff=0.1, order=5):
    fs = fs1*len(data.iloc[:,1])
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data.iloc[:, 1])
    time_vec = data.iloc[:, 0]
    smooth_df = pd.DataFrame({'TIME': time_vec, 'FAM': y})
    return smooth_df

def plot_data(dfinit, df0, df1, df2, filename, th):
    max1, t1 = find_1st(df1)
    max2, t2 = find_2nd(df2)
    th_time, th_ind = find_threshold(df0, th)
    actual_time = df0.iloc[th_time, 0]
    print(df0.iloc[th_time, 1], th_time)

    fig, ax = plt.subplots(2,2)
    ax[0, 0].plot(dfinit.iloc[:, 0], dfinit.iloc[:, 1])
    ax[0, 0].title.set_text('Native')
    ax[0,1].plot(df0.iloc[:, 0], df0.iloc[:, 1])
    ax[0,1].plot(t1*np.ones(10), np.linspace(0, df0.iloc[-1, 1],10), color='g')
    ax[0,1].plot(t2*np.ones(10), np.linspace(0, df0.iloc[-1, 1],10), color='r')
    ax[0, 1].plot(np.linspace(0, df0.iloc[-1, 0],10), th*np.ones(10))
    ax[0, 1].scatter(actual_time, df0.iloc[th_time, 1], color='b')
    ax[0,1].title.set_text('Filtered')
    ax[1,0].plot(df1.iloc[:, 0], df1.iloc[:, 1])
    ax[1,0].scatter(t1, max1, color='g')
    ax[1,0].title.set_text('1st derivative')
    ax[1,1].plot(df2.iloc[:, 0], df2.iloc[:, 1])
    ax[1,1].scatter(t2, max2, color='r')
    ax[1,1].title.set_text('2nd derivative')
    fig.tight_layout()
    fig.savefig(filename + '.png')
    plt.clf()
    return

def find_1st(df):
    list = df.iloc[:,1]
    max = 0
    time = 0
    #maximum needs to be between minute 2 and minute 10 of run
    time_wind_low = 0.1*len(list)
    time_wind_up = 1*len(list)

    for num, val in enumerate(list):
        if val > max and int(time_wind_low) < num < int(time_wind_up):
            max = val
            time = df.iloc[num, 0]

    return max, time

def find_2nd(df):
    list = df.iloc[:,1]
    max = 0
    time = 0
    #maximum needs to be between minute 2 and minute 10 of run
    time_wind_low = 0.1*len(list)
    time_wind_up = 0.5*len(list)

    for num, val in enumerate(list):
        if val > max and int(time_wind_low) < num < int(time_wind_up):
            max = val
            time = df.iloc[num, 0]

    return max, time

def find_threshold(df, th):
    indices = [i for i,v in enumerate(df.iloc[:,1]) if v > th]
    ind_th = 0
    try:
        th_time = min(indices)
    except:
        th_time = 0
        ind_th = 0
    return th_time, ind_th

if __name__ == "__main__":
    save_dir = "Results_Storage/Plots"
    os.makedirs(save_dir, exist_ok=True)

    filt_param = [4, 10]
    fs = filt_param[1]
    resultdat = []
    filelist = []

    #prediction parameters
    th = 0.03
    cutoff_T1 = 60
    cutoff_T2_low = 60
    cutoff_T2_high = 600

    counter = 1

    for dirpath, dirnames, filenames in os.walk("./INMI_CSV", topdown=True):
        for f in filenames:
            if f.endswith('_pcameasurements.dat'):
                filelist.append(os.path.join(dirpath, f))
                counter = counter + 1

    print("Number of files: ", counter)

    for i, file in enumerate(filelist[0:250]):
        splitpath = file.split("\\")
        exp_dat = splitpath[1].split("_")
        filename_img = splitpath[-1]

        print(i)
        data = encoding_read(file)
        filt_factor1 = filt_param[0] * round(0.01 * len(data))
        filt_factor2 = filt_param[1] * round(0.01 * len(data))

        smooth_data = smooth_filt(data, filt_factor1)
        grad = [gradient1(smooth_data), gradient1(smooth_filt(gradient1(smooth_data), filt_factor1))]
        grad_smooth = [smooth_filt(grad[0], filt_factor2), smooth_filt(grad[1], filt_factor2)]

        plot_data(data, smooth_data, grad_smooth[0], grad_smooth[1], save_dir + "\\" + filename_img, th)

        maxval1, time1 = find_1st(grad_smooth[0])
        maxval2, time2 = find_2nd(grad_smooth[1])

        Tth, ind_th = find_threshold(smooth_data, th)

        # print(maxval, maxind)
        if float(max(data.iloc[:, 1])) > th and time1 > cutoff_T1 and time2 > cutoff_T2_low and time2 < cutoff_T2_high:
            prediction = "yes"
        else:
            prediction = "no"
        results = [exp_dat[-3], exp_dat[-2], exp_dat[-1], float(max(data.iloc[:, 1])), maxval1, time1, maxval2, time2, Tth,
                prediction, file]
        resultdat.append(results)

# plot_data(data, smooth_data, grad_filt[0], grad_filt[1])
restab = pd.DataFrame(resultdat, columns=['Date', 'NanoID', 'SampleID', 'MAX RFU', '1der MAX', 'T1', '2der MAX', 'T2', 'Tth', 'prediction', 'filename'])
print('DATAframe generated')

restab.to_csv(save_dir + "/" + "240606_Pproc_INMI.csv", encoding='utf-8')
print('Files stored to: ', save_dir)
print(restab.head())
