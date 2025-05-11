#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:57:18 2024

@author: fmuir
"""


import os
import pickle
import datetime as dt
from datetime import datetime
import time
import numpy as np
import random

import pandas as pd
import geopandas as gpd
pd.options.mode.chained_assignment = None # suppress pandas warning about setting a value on a copy of a slice
from scipy.interpolate import interp1d, PchipInterpolator
from scipy.stats import circmean

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering
# from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.metrics import silhouette_score, root_mean_squared_error
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import train_test_split
# from sklearn.utils.class_weight import compute_class_weight
# from imblearn.over_sampling import SMOTE
from tensorflow.keras import backend as K
from tensorflow.keras.models import Sequential
from tensorflow.keras import Input
from tensorflow.keras.layers import LSTM, Dense, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.optimizers.schedules import ExponentialDecay
from imblearn.over_sampling import SMOTE
from tensorflow.keras.callbacks import EarlyStopping, TensorBoard
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator
from tensorboard.plugins.hparams import api as hp

import optuna
import shap

# Only use tensorflow in CPU mode
import tensorflow as tf
tf.config.set_visible_devices([],'GPU')


# ----------------------------------------------------------------------------------------
#%% DATA PREPPING FUNCTIONS ###

def LoadIntersections(filepath, sitename):
    """
    Load in transect intersection dataframes stored as pickle files. Generated 
    using Transects.GetIntersections(), Transects.SaveIntersections(), 
    Transects.GetBeachWidth(), Transects.SaveWaterIntersections(),
    Transects.TZIntersect(), Transects.SlopeIntersect(), 
    Transects.WavesIntersect().
    FM Jul 2024

    Parameters
    ----------
    filepath : str
        Path to 'Data' directory for chosen site.
    sitename : str
        Name of site chosen.

    Returns
    -------
    TransectInterGDF : GeoDataFrame
        GeoDataFrame of cross-shore transects, intersected with vegetation edge lines.
    TransectInterGDFWater : GeoDataFrame
        GeoDataFrame of cross-shore transects, intersected with waterlines.
    TransectInterGDFTopo : GeoDataFrame
        GeoDataFrame of cross-shore transects, intersected with slope raster and vegetation transition zones.
    TransectInterGDFWave : GeoDataFrame
        GeoDataFrame of cross-shore transects, intersected with Copernicus hindcast wave data.

    """
    print('Loading transect intersections...')
    
    with open(os.path.join
              (filepath , sitename, 'intersections', sitename + '_transect_intersects.pkl'), 'rb') as f:
        TransectInterGDF = pickle.load(f)
        
    with open(os.path.join
              (filepath , sitename, 'intersections', sitename + '_transect_water_intersects.pkl'), 'rb') as f:
        TransectInterGDFWater = pickle.load(f)

    with open(os.path.join
              (filepath , sitename, 'intersections', sitename + '_transect_topo_intersects.pkl'), 'rb') as f:
        TransectInterGDFTopo = pickle.load(f)

    with open(os.path.join
              (filepath , sitename, 'intersections', sitename + '_transect_wave_intersects.pkl'), 'rb') as f:
        TransectInterGDFWave = pickle.load(f)
        
    return TransectInterGDF, TransectInterGDFWater, TransectInterGDFTopo, TransectInterGDFWave
        

def CompileTransectData(TransectInterGDF, TransectInterGDFWater, TransectInterGDFTopo, TransectInterGDFWave):
    """
    Merge together transect geodataframes produced from COASTGUARD.VedgeSat and CoastSat. Each transect holds 
    timeseries of a range of satellite-derived metrics.
    FM Aug 2024

    Parameters
    ----------
    TransectInterGDF : GeoDataFrame
        DataFrame of cross-shore transects intersected with timeseries of veg edges.
    TransectInterGDFWater : GeoDataFrame
        DataFrame of cross-shore transects intersected with timeseries of waterlines.
    TransectInterGDFTopo : GeoDataFrame
        DataFrame of cross-shore transects intersected with timeseries of slopes at the veg edge.
    TransectInterGDFWave : GeoDataFrame
        DataFrame of cross-shore transects intersected with timeseries of wave conditions.

    Returns
    -------
    CoastalDF : DataFrame
        DataFrame processed.

    """
    print('Merging transect-based data...')
    # Merge veg edge intersection data with waterline intersection data
    CoastalDF = pd.merge(TransectInterGDF[['TransectID',
                                           'dates','times','distances']], 
                         TransectInterGDFWater[['TransectID',
                                                'wldates','wltimes','wlcorrdist', 
                                                'tideelev','beachwidth',
                                                'tidedatesFD','tideelevFD', 'tideelevMx']],
                         how='inner', on='TransectID')
    # Merge combined dataframe with wave info
    # TransectInterGDFWave[['TransectID','WaveHs', 'WaveDir', 'WaveTp', 'WaveDiffus']]
    # CoastalDF = pd.merge(CoastalDF, 
    #                      TransectInterGDFWave[['TransectID','WaveDates','WaveHs', 'WaveDir', 'WaveTp', 'WaveAlpha', 'Runups','Iribarrens']],
    #                      how='inner', on='TransectID')
    CoastalDF = pd.merge(CoastalDF, 
                         TransectInterGDFWave[['TransectID',
                                               'WaveDates','WaveDatesFD',
                                               'WaveHsFD', 'WaveDirFD', 'WaveTpFD', 'WaveAlphaFD', 
                                               'Runups', 'Iribarren']],
                         how='inner', on='TransectID')
    
    CoastalDF.rename(columns={'distances':'VE', 'wlcorrdist':'WL'}, inplace=True)
    
    print('Converting to datetimes...')
    veDTs = []
    for Tr in range(len(CoastalDF)):
        veDTs_Tr = []
        for i in range(len(CoastalDF['dates'].iloc[Tr])):
            veDTs_Tr.append(dt.datetime.strptime(CoastalDF['dates'].iloc[Tr][i]+' '+CoastalDF['times'].iloc[Tr][i], '%Y-%m-%d %H:%M:%S.%f'))
        veDTs.append(veDTs_Tr)
    CoastalDF['veDTs'] = veDTs
    
    wlDTs = []
    for Tr in range(len(CoastalDF)):
        wlDTs_Tr = []
        for i in range(len(CoastalDF['wldates'].iloc[Tr])):
            wlDTs_Tr.append(dt.datetime.strptime(CoastalDF['wldates'].iloc[Tr][i]+' '+CoastalDF['wltimes'].iloc[Tr][i], '%Y-%m-%d %H:%M:%S.%f'))
        wlDTs.append(wlDTs_Tr)
    CoastalDF['wlDTs'] = wlDTs
    
    
    return CoastalDF


def InterpWL(CoastalDF, Tr):
    """
    Interpolate over waterline associated timeseries so that dates 
    match vegetation associated ones.
    NOT USED
    FM Aug 2024

    Parameters
    ----------
    CoastalDF : DataFrame
        DataFrame of cross-shore transects (rows) and intersected coastal 
        timeseries/metrics (columns).
    Tr : int
        Transect ID of choice.

    Returns
    -------
    TransectDF : DataFrame
        Subset row matching the requested transect ID (Tr), with interpolated
        values for 'wlcorrdist', 'waterelev' and 'beachwidth'.

    """
    TransectDF = CoastalDF.iloc[[Tr],:] # single-row dataframe
    # TransectDF = TransectDF.transpose()

    # Interpolate over waterline associated variables to match dates with veg edge dates
    wl_numdates = pd.to_datetime(TransectDF['wldates'][Tr]).values.astype(np.int64)
    ve_numdates = pd.to_datetime(TransectDF['dates'][Tr]).values.astype(np.int64)
    wl_interp_f = interp1d(wl_numdates, TransectDF['wlcorrdist'][Tr], kind='linear', fill_value='extrapolate')
    wl_interp = wl_interp_f(ve_numdates).tolist()
    welev_interp_f = interp1d(wl_numdates, TransectDF['tideelev'][Tr], kind='linear', fill_value='extrapolate')
    welev_interp = welev_interp_f(ve_numdates).tolist()
    TransectDF['wlcorrdist'] = [wl_interp]
    TransectDF['tideelev'] = [welev_interp]
    # Recalculate beachwidth
    beachwidth = [abs(wl_interp[i] - TransectDF['distances'][Tr][i]) for i in range(len(wl_interp))]
    TransectDF['beachwidth'] = [beachwidth]
    
    TransectDF.drop(columns=['wldates'], inplace=True)
    
    # Transpose to get columns of variables and rows of timesteps
    TransectDF = pd.DataFrame({col: pd.Series(val.iloc[0]) for col,val in TransectDF.items()})
    
    return TransectDF


def InterpWLWaves(CoastalDF, Tr):
    """
    Interpolate over waterline and wave associated timeseries so that dates 
    match vegetation associated ones.
    NOT USED
    FM Aug 2024

    Parameters
    ----------
    CoastalDF : DataFrame
        DataFrame of cross-shore transects (rows) and intersected coastal 
        timeseries/metrics (columns).
    Tr : int
        Transect ID of choice.

    Returns
    -------
    TransectDF : DataFrame
        Subset row matching the requested transect ID (Tr), with interpolated
        values for 'wlcorrdist', 'waterelev' and 'beachwidth'.

    """
    TransectDF = CoastalDF.iloc[[Tr],:] # single-row dataframe
    # TransectDF = TransectDF.transpose()

    # Interpolate over waterline and wave associated variables
    wl_numdates = pd.to_datetime(TransectDF['wlDTs'][Tr]).values.astype(np.int64)
    ve_numdates = pd.to_datetime(TransectDF['veDTs'][Tr]).values.astype(np.int64)
    wv_numdates = pd.to_datetime(TransectDF['WaveDates'][Tr]).values.astype(np.int64)
    # Match dates with veg edge dates and append back to TransectDF
    for wlcol in ['wlcorrdist', 'tideelev','beachwidth']:
        wl_interp_f = interp1d(wl_numdates, TransectDF[wlcol][Tr], kind='linear', fill_value='extrapolate')
        wl_interp = wl_interp_f(ve_numdates).tolist()
        TransectDF[wlcol] = [wl_interp]
    for wvcol in ['WaveHs','WaveDir','WaveTp','Runups', 'Iribarrens']:
        wv_interp_f = interp1d(wv_numdates, TransectDF[wvcol][Tr], kind='linear', fill_value='extrapolate')
        wv_interp = wv_interp_f(ve_numdates).tolist()
        TransectDF[wvcol] = [wv_interp]
    
    # Recalculate beachwidth as values will now be mismatched
    beachwidth = [abs(wl_interp[i] - TransectDF['distances'][Tr][i]) for i in range(len(wl_interp))]
    TransectDF['beachwidth'] = [beachwidth]
    
    TransectDF.drop(columns=['WaveDates','wldates','wltimes', 'wlDTs','dates','times'], inplace=True)
    
    # Transpose to get columns of variables and rows of timesteps
    TransectDF = pd.DataFrame({col: pd.Series(val.iloc[0]) for col,val in TransectDF.items()})
    
    # Reset index for timeseries
    TransectDF.index = TransectDF['veDTs']
    TransectDF = TransectDF.drop(columns=['TransectID', 'veDTs'])

    return TransectDF


def InterpVEWLWv(CoastalDF, Tr, IntpKind='linear'):
    """
    Interpolate over waterline and vegetation associated timeseries so that dates 
    match wave associated (full timeseries) ones.
    FM Nov 2024

    Parameters
    ----------
    CoastalDF : DataFrame
        DataFrame of cross-shore transects (rows) and intersected coastal 
        timeseries/metrics (columns).
    Tr : int
        Transect ID of choice.

    Returns
    -------
    TransectDF : DataFrame
        Subset row matching the requested transect ID (Tr), with interpolated
        values for 'wlcorrdist', 'waterelev' and 'beachwidth'.

    """
    TransectDF = CoastalDF.iloc[[Tr],:] # single-row dataframe
    # TransectDF = TransectDF.transpose()

    # Convert all circular metrics to radians
    for WaveProp in ['WaveDirFD','WaveAlphaFD']:
        TransectDF[WaveProp] = [np.deg2rad(WaveDir) for WaveDir in TransectDF[WaveProp]]

    # Sort and remove duplicate dates
    def Preprocess(dates, values):
        df = pd.DataFrame({'dates': dates, 'values': values})
        # Group by unique dates and take the mean of duplicate values
        df_grouped = df.groupby('dates', as_index=False).mean()
        # If any nans still exist, fill them with PCHIP
        if np.isnan(df_grouped['values']).sum() > 0:
            df_grouped = df_grouped.interpolate(kind='pchip')
        # Extract clean dates and values
        unique_numdates = df_grouped['dates'].values
        unique_vals = df_grouped['values'].values
        
        return unique_numdates, unique_vals
    
    # Calculate circular exponential weighted average
    def CircEWA(DFrad):
        # Compute sine and cosine
        DFsin = np.sin(DFrad)
        DFcos = np.cos(DFrad)
        
        # Apply EWA to sine and cosine components
        DFEWAsin = DFsin.ewm(span=10, adjust=False).mean()
        DFEWAcos = DFcos.ewm(span=10, adjust=False).mean()
        
        return np.arctan2(DFEWAsin, DFEWAcos)


    # Interpolate over waterline and wave associated variables
    wl_numdates = pd.to_datetime(TransectDF['wlDTs'][Tr]).values.astype(np.int64)
    ve_numdates = pd.to_datetime(TransectDF['veDTs'][Tr]).values.astype(np.int64)
    td_numdates = pd.to_datetime(TransectDF['tidedatesFD'][Tr]).values.astype(np.int64)
    wvsat_numdates = pd.to_datetime(TransectDF['WaveDates'][Tr]).values.astype(np.int64)
    wv_numdates = pd.to_datetime(TransectDF['WaveDatesFD'][Tr]).values.astype(np.int64)
    # Adjacent VE and WL positions (+1 is to the upcoast of current transect, -1 is to downcoast)
    wl_numdates_up = pd.to_datetime(CoastalDF.iloc[[Tr+1],:]['wlDTs'][Tr+1]).values.astype(np.int64)
    ve_numdates_up = pd.to_datetime(CoastalDF.iloc[[Tr+1],:]['veDTs'][Tr+1]).values.astype(np.int64)
    wl_numdates_down = pd.to_datetime(CoastalDF.iloc[[Tr-1],:]['wlDTs'][Tr-1]).values.astype(np.int64)
    ve_numdates_down = pd.to_datetime(CoastalDF.iloc[[Tr-1],:]['veDTs'][Tr-1]).values.astype(np.int64)
    
    for wvcol in ['WaveDirFD', 'WaveAlphaFD', 'WaveHsFD', 'WaveTpFD']:
        if IntpKind == 'pchip':
            wv_numdates_clean, wvcol_clean = Preprocess(wv_numdates, TransectDF[wvcol][Tr])
            wv_interp_f = PchipInterpolator(wv_numdates_clean, wvcol_clean)
        else:
            wv_interp_f = interp1d(wv_numdates, TransectDF[wvcol][Tr], kind=IntpKind, fill_value='extrapolate')
        wv_interp = wv_interp_f(wv_numdates).tolist()
        TransectDF[wvcol] = [wv_interp]
        
    # Match dates with veg edge dates and append back to TransectDF
    for wlcol in ['WL', 'tideelev','beachwidth']:
        if IntpKind == 'pchip':
            wl_numdates_clean, wlcol_clean = Preprocess(wl_numdates, TransectDF[wlcol][Tr])
            wl_interp_f = PchipInterpolator(wl_numdates_clean, wlcol_clean)
        else:
            wl_interp_f = interp1d(wl_numdates, TransectDF[wlcol][Tr], kind=IntpKind, fill_value='extrapolate')
        wl_interp = wl_interp_f(wv_numdates).tolist()
        TransectDF[wlcol] = [wl_interp]
    for wvsatcol in ['Runups','Iribarren']:
        if IntpKind == 'pchip':
            wvsat_numdates_clean, wvsatcol_clean = Preprocess(wvsat_numdates, TransectDF[wvsatcol][Tr])
            wvsat_interp_f = PchipInterpolator(wvsat_numdates_clean, wvsatcol_clean)
        else:
            wvsat_interp_f = interp1d(wvsat_numdates, TransectDF[wvsatcol][Tr], kind=IntpKind, fill_value='extrapolate')
        wvsat_interp = wvsat_interp_f(wv_numdates).tolist()
        TransectDF[wvsatcol] = [wvsat_interp]
    for tdcol in ['tideelevFD','tideelevMx']:
        if IntpKind == 'pchip':
            td_numdates_clean, tdcol_clean = Preprocess(td_numdates, TransectDF[tdcol][Tr])
            td_interp_f = PchipInterpolator(td_numdates_clean, tdcol_clean)
        else:
            td_interp_f = interp1d(td_numdates, TransectDF[tdcol][Tr], kind=IntpKind, fill_value='extrapolate')
        td_interp = td_interp_f(wv_numdates).tolist()
        TransectDF[tdcol] = [td_interp]
    for vecol in ['VE']:#,'TZwidth']:
        if IntpKind == 'pchip':
            ve_numdates_clean, vecol_clean = Preprocess(ve_numdates, TransectDF[vecol][Tr])
            ve_interp_f = PchipInterpolator(ve_numdates_clean, vecol_clean)
        else:
            ve_interp_f = interp1d(ve_numdates, TransectDF[vecol][Tr], kind=IntpKind, fill_value='extrapolate')
        ve_interp = ve_interp_f(wv_numdates).tolist()
        TransectDF[vecol] = [ve_interp]
    
    # Adjacent VE and WL interpolation
    for wlcol in ['WL']:
        if IntpKind == 'pchip':
            wl_numdates_clean, wlcol_clean = Preprocess(wl_numdates_up, CoastalDF.iloc[[Tr+1],:][wlcol][Tr+1])
            wl_interp_f = PchipInterpolator(wl_numdates_clean, wlcol_clean)
        else:
            wl_interp_f = interp1d(wl_numdates_up, CoastalDF.iloc[[Tr+1],:][wlcol][Tr+1], kind=IntpKind, fill_value='extrapolate')
        wl_interp = wl_interp_f(wv_numdates).tolist()
        TransectDF[wlcol+'_u'] = [wl_interp]
    for vecol in ['VE']:#,'TZwidth']:
        if IntpKind == 'pchip':
            ve_numdates_clean, vecol_clean = Preprocess(ve_numdates_up, CoastalDF.iloc[[Tr+1],:][vecol][Tr+1])
            ve_interp_f = PchipInterpolator(ve_numdates_clean, vecol_clean)
        else:
            ve_interp_f = interp1d(ve_numdates_up, CoastalDF.iloc[[Tr+1],:][vecol][Tr+1], kind=IntpKind, fill_value='extrapolate')
        ve_interp = ve_interp_f(wv_numdates).tolist()
        TransectDF[vecol+'_u'] = [ve_interp]
    for wlcol in ['WL']:
        if IntpKind == 'pchip':
            wl_numdates_clean, wlcol_clean = Preprocess(wl_numdates_down, CoastalDF.iloc[[Tr-1],:][wlcol][Tr-1])
            wl_interp_f = PchipInterpolator(wl_numdates_clean, wlcol_clean)
        else:
            wl_interp_f = interp1d(wl_numdates_down, CoastalDF.iloc[[Tr-1],:][wlcol][Tr-1], kind=IntpKind, fill_value='extrapolate')
        wl_interp = wl_interp_f(wv_numdates).tolist()
        TransectDF[wlcol+'_d'] = [wl_interp]
    for vecol in ['VE']:#,'TZwidth']:
        if IntpKind == 'pchip':
            ve_numdates_clean, vecol_clean = Preprocess(ve_numdates_down, CoastalDF.iloc[[Tr-1],:][vecol][Tr-1])
            ve_interp_f = PchipInterpolator(ve_numdates_clean, vecol_clean)
        else:
            ve_interp_f = interp1d(ve_numdates_down, CoastalDF.iloc[[Tr-1],:][vecol][Tr-1], kind=IntpKind, fill_value='extrapolate')
        ve_interp = ve_interp_f(wv_numdates).tolist()
        TransectDF[vecol+'_d'] = [ve_interp]
    
    # Recalculate beachwidth as values will now be mismatched
    beachwidth = [abs(wl_interp[i] - TransectDF['VE'][Tr][i]) for i in range(len(wl_interp))]
    TransectDF['beachwidth'] = [beachwidth]
    
    TransectDF.drop(columns=['veDTs','wldates','wltimes', 'wlDTs','dates','times'], inplace=True)
    
    # Transpose to get columns of variables and rows of timesteps
    TransectDF = pd.DataFrame({col: pd.Series(val.iloc[0]) for col,val in TransectDF.items()})
    
    # Reset index for timeseries
    TransectDF.index = TransectDF['WaveDatesFD']
    TransectDF = TransectDF.drop(columns=['TransectID', 'WaveDates','WaveDatesFD', 'tidedatesFD'])

    # If any index rows (i.e. daily wave dates) are empty, remove them
    TransectDF = TransectDF[~pd.isnull(TransectDF.index)]
    
    # Convert circular variables to components (to conserve topology and avoid jumps)
    TransectDF['WaveDirsin'] = np.sin(TransectDF['WaveDirFD'])
    TransectDF['WaveDircos'] = np.cos(TransectDF['WaveDirFD'])
    
    # Calculate moving window averages (MA) and exponentially weighted averages (EW) for each daily wave condition
    for wvcol in ['WaveHsFD', 'WaveTpFD']:
        TransectDF[wvcol[:-2]+'MA'] = TransectDF[wvcol].rolling(window=10, min_periods=1).mean()
        TransectDF[wvcol[:-2]+'EW'] = TransectDF[wvcol].ewm(span=10, adjust=False).mean()
    for wvcol in ['WaveDirFD', 'WaveAlphaFD']:

        TransectDF[wvcol[:-2]+'MA'] = TransectDF[wvcol].rolling(window=10, min_periods=1).apply(
                                    lambda x: circmean(x, high=360, low=0), raw=True)
        TransectDF[wvcol[:-2]+'EW'] = CircEWA(TransectDF[wvcol])
        TransectDF[wvcol[:-2]+'sinEW'] = np.sin(TransectDF[wvcol[:-2]+'EW'])
        TransectDF[wvcol[:-2]+'cosEW'] = np.cos(TransectDF[wvcol[:-2]+'EW'])
    
    # Time-lagging neighbour values by 10 days (and fill gap at start with mean value)
    for col in ['WL_u', 'VE_u', 'WL_d', 'VE_d']:
        TransectDF[col+'-10'] = TransectDF[col].shift(periods=10, axis=0, fill_value=TransectDF[col].mean())
    
    return TransectDF


def DailyInterp(TransectDF):
    """
    Reindex and interpolate over timeseries to convert it to daily
    FM Nov 2024

    Parameters
    ----------
    VarDF : DataFrame
        Dataframe of per-transect coastal metrics/variables in (irregular) timeseries.

    Returns
    -------
    VarDFDay : DataFrame
        Dataframe of per-transect coastal metrics/variables in daily timesteps.
    """
    
    # Fill nans factoring in timesteps for interpolation
    TransectDF.replace([np.inf, -np.inf], np.nan, inplace=True)
    VarDF = TransectDF.interpolate(method='time', axis=0)
    # Strip time (just set to date timestamps), and take mean if any timesteps duplicated
    VarDF.index = VarDF.index.normalize()
    VarDF = VarDF.groupby(VarDF.index).mean()
    
    DailyInd = pd.date_range(start=VarDF.index.min(), end=VarDF.index.max(), freq='D')
    VarDFDay = VarDF.reindex(DailyInd)
    
    for col in VarDFDay.columns:
        if col == 'WaveDir' or col == 'WaveAlpha':
            # VarDFDay[col] = VarDFDay[col].interpolate(method='pchip')
            # Separate into components to allow interp over 0deg threshold
            # but need to preserve the offset between WaveDir and WaveAlpha
            VarDFsin = np.sin(VarDFDay[col])
            VarDFsininterp = VarDFsin.interpolate(method='pchip')
            VarDFcos = np.cos(VarDFDay[col])
            VarDFcosinterp = VarDFcos.interpolate(method='pchip')
            VarDFDay[col] = np.arctan2(VarDFsininterp, VarDFcosinterp) % (2 * np.pi) # keep values 0-2pi
        else:
            VarDFDay[col] = VarDFDay[col].interpolate(method='spline',order=1)
    
    return VarDFDay


def CreateSequences(X, y=None, time_steps=1, Pad=False):
    '''
    Function to create sequences (important for timeseries data where data point
    is temporally dependent on the one that came before it). Data sequences are
    needed for training RNNs, where temporal patterns are learned.
    FM June 2024

    Parameters
    ----------
    X : array
        Training data as array of feature vectors.
    y : array
        Training classes as array of binary labels.
    time_steps : int, optional
        Number of time steps over which to generate sequences. The default is 1.
    Pad : bool, optional
        Flag for whether to apply a padding of zeros to the start of the timeseries
        if there isn't enough data for a full sequence. The default is False.

    Returns
    -------
    array, array
        Numpy arrays of sequenced data

    '''
    # Initialise training and target feature lists
    Xs = []
    ys = []
    Ind = []
    # X_array = X.to_numpy(dtype='float32')  
    # y_array = y.to_numpy(dtype='float32') if y is not None else None
    # if Pad:
    if len(X) < time_steps:
        pad_length = time_steps - len(X)+1
        pad_array = np.zeros((pad_length, X.shape[1]), dtype='float32')
        Xarr = np.vstack((X.to_numpy(dtype='float32'), pad_array))
        if y is not None:
            pad_array = np.zeros((pad_length, y.shape[1]), dtype='float32')
            yarr = np.vstack((y.to_numpy(dtype='float32'), pad_array))
    else:
        if y is not None:
            Xarr, yarr = X.to_numpy(dtype='float32'), y.to_numpy(dtype='float32')
        else:
            Xarr, yarr = X.to_numpy(dtype='float32'), None
    if len(Xarr) >= time_steps:  # Check if there's enough data
        for i in range(len(Xarr) - time_steps):
            # Slice feature set into sequences using moving window of size = number of timesteps
            Xs.append(Xarr[i:(i + time_steps)])
            Ind.append(X.index[min(i + time_steps, len(X) - 1)])  # Preserve index
            if y is not None and i + time_steps < len(y):
                ys.append(yarr[i + time_steps])
            else:
                ys.append(np.nan)
        return np.array(Xs), np.array(ys), np.array(Ind)
    else:
        # Not enough data to create a sequence
        print(f"Not enough data to create sequences with time_steps={time_steps}")
        return np.array([]), np.array([]), np.array([])


# ----------------------------------------------------------------------------------------
#%% MODEL INFRASTRUCTURE FUNCTIONS ###


def CostSensitiveLoss(falsepos_cost, falseneg_cost, binary_thresh):
    """
    Create a cost-sensitive loss function to implement within an NN model.compile() step.
    NN model is penalised based on cost matrix, which helps avoid mispredictions
    with disproportionate importance (i.e. a missed storm warning is more serious
    than a false alarm).
    FM June 2024

    Parameters
    ----------
    falsepos_cost : int
        Proportional weight towards false positive classification.
    falseneg_cost : int
        Proportional weight towards false negative classification.
    binary_thresh : float
        Value between 0 and 1 representing .

    Returns
    -------
    loss : function
        Calls the loss function when it is set within model.compile(loss=LossFn).

    """
    def loss(y_true, y_pred):
        # Flatten the arrays
        y_true = K.flatten(y_true)
        y_pred = K.flatten(y_pred)
        
        # Convert predictions to binary class predictions
        y_pred_classes = K.cast(K.greater(y_pred, binary_thresh), tf.float32)
        
        # Calculate cost
        falsepos = K.sum(y_true * (1 - y_pred_classes) * falseneg_cost)
        falseneg = K.sum((1 - y_true) * y_pred_classes * falsepos_cost)
        
        bce = K.binary_crossentropy(y_true, y_pred)
        
        return bce + falsepos + falseneg
    
    return loss


@tf.keras.utils.register_keras_serializable()
def ShoreshopLoss(y_true, y_pred):
    """
    Apply a loss function to the shoreline prediction model using the Shoreshop
    guidelines function, which is a mix of the RMSE (normalised), Pearson 
    Correlation and Standard Deviation (normalised). More info at:
    https://github.com/ShoreShop/ShoreModel_Benchmark/tree/main
    
    Note: tf math functions are used because they return tensors, not arrays or
    scalars. This is so GPU acceleration and backpropagation can take place.
    FM Mar 2025
    
    Parameters
    ----------
    y_true : float, array
        Actual values.
    y_pred : float, array
        Predicted values to measure against actual values.

    Returns
    -------
    loss : float
        Measure of error based on actual vs predicted.

    """
    # Calculate RMSE
    rmse_pred = tf.sqrt(tf.reduce_mean(tf.square(y_true - y_pred)))
    
    # Calculate Standard Deviation
    std_targ = tf.math.reduce_std(y_true)
    std_pred = tf.math.reduce_std(y_pred)
    
    # Calculate normalised versions of RMSE and StDv
    rmse_norm = rmse_pred / (std_targ + 1e-8)  # Avoid division by zero
    std_norm = std_pred / (std_targ + 1e-8)
    
    # Calculate Correlation Coefficient
    mean_true = tf.reduce_mean(y_true)
    mean_pred = tf.reduce_mean(y_pred)
    
    num = tf.reduce_sum((y_true - mean_true) * (y_pred - mean_pred))
    den = tf.sqrt(tf.reduce_sum(tf.square(y_true - mean_true)) * tf.reduce_sum(tf.square(y_pred - mean_pred)) + 1e-8)    
    corr = num / den  # Pearson correlation coefficient
    
    # Calculate final loss 
    loss = tf.sqrt(tf.square(0 - rmse_norm) + tf.square(1 - corr) + tf.square(1 - std_norm))
    
    return loss


def PrepData(TransectDF, MLabels, ValidSizes, TSteps, TrainFeatCols, TargFeatCols, TrainTestPortion=None, UseSMOTE=False):
    """
    Prepare features (X) and labels (y) for feeding into a NN for timeseries prediction.
    Timeseries is expanded and interpolated to get daily timesteps, then scaled across
    all variables, then split into model inputs and outputs before sequencing these
    
    FM Sept 2024
    
    Parameters
    ----------
    TransectDF : DataFrame
        Dataframe of per-transect coastal metrics/variables in timeseries.
    MLabels : list
        List of unique model run names.
    TestSizes : list, hyperparameter
        List of different proportions of data to separate out from training to validate model.
    TSteps : list, hyperparameter
        Number of timesteps (1 step = 1 day) to sequence training data over.
    UseSMOTE : bool, optional
        Flag for using SMOTE to oversample imbalanced data. The default is False.

    Returns
    -------
    PredDict : dict
        Dictionary to store all the NN model metadata.
    VarDFDay_scaled : DataFrame
        Scaled DataFrame of past data interpolated to daily timesteps (with temporal index), 
        for training and validation.
    VarDFDayTest_scaled : DataFrame
        Scaled DataFrame of past data interpolated to daily timesteps (with temporal index), 
        for testing (unseen by model).

    """
    
    # Interpolate over data gaps to get regular (daily) measurements
    # VarDFDay = DailyInterp(TransectDF)
    
    
    
    # Scale vectors to normalise themTrainTestPortion is not None and 
    Scalings = {}
    TransectDF_scaled = TransectDF.copy()
    for col in TransectDF.columns:
        Scalings[col] = StandardScaler()
        TransectDF_scaled[col] = Scalings[col].fit_transform(TransectDF[[col]])
      
    if TrainTestPortion is not None and type(TrainTestPortion) == float:
        # if train-test proportion is a percentage
        VarDFDay_scaled = TransectDF_scaled.iloc[:int(len(TransectDF_scaled)*TrainTestPortion)]
        VarDFDayTest_scaled = TransectDF_scaled.iloc[int(len(TransectDF_scaled)*TrainTestPortion):]
    elif TrainTestPortion is not None and type(TrainTestPortion) == datetime:
        # if the train-test proportion is a date to split between
        VarDFDay_scaled = TransectDF_scaled.loc[:TrainTestPortion]
        VarDFDayTest_scaled = TransectDF_scaled.loc[TrainTestPortion:]
    else:
        # just do last 20%
        VarDFDay_scaled = TransectDF_scaled.iloc[:int(len(TransectDF_scaled)*0.8)]
        VarDFDayTest_scaled = TransectDF_scaled.iloc[int(len(TransectDF_scaled)*0.8):]
    
    # Define prediction dictionary for multiple runs/hyperparameterisation
    PredDict = {'mlabel':MLabels,   # name of the model run
                'model':[],         # compiled model (Sequential object)
                'history':[],       # training history to be written to later
                'loss':[],          # final loss value of run
                'accuracy':[],      # final accuracy value of run
                'train_time':[],    # time taken to train
                'seqlen':[],        # length of temporal sequence in timesteps to break data up into
                'trainfeats':[],    # column names of training features
                'targfeats':[],     # column names of target features
                'validsize':[],     # portion (decimal) of training data to set aside as validation 
                'X_train':[],       # training features/cross-shore values (scaled and filled to daily)
                'y_train':[],       # training target/cross-shore VE and WL (scaled and filled to daily)
                'X_val':[],         # validation features/cross-shore values
                'y_val':[],         # validation target/cross-shore VE and WL
                'scalings':[],      # scaling used on each feature (to transform them back to real values)
                'epochN':[],        # number of times the full training set is passed through the model
                'batchS':[],        # number of samples to use in one iteration of training
                'denselayers':[],   # number of dense layers in model construction
                'dropoutRt':[],     # percentage of nodes to randomly drop in training (avoids overfitting)
                'learnRt':[],       # size of steps to adjust parameters by on each iteration
                'hiddenLscale':[]}  # Scaling relationship for number of hidden LSTM layers 
    
    for MLabel, ValidSize, TStep, TrainFeatCol, TargFeatCol in zip(PredDict['mlabel'], ValidSizes, TSteps, TrainFeatCols, TargFeatCols):
        
        # Separate into training features (what will be learned from) and target features (what will be predicted)
        TrainFeat = VarDFDay_scaled[TrainFeatCol]
        # TrainFeat = VarDFDay_scaled[['WaveDir', 'Runups', 'Iribarren']]
        TargFeat = VarDFDay_scaled[TargFeatCol] # vegetation edge and waterline positions
        
        # Add scaling relationships and sequence lengths to dict to convert back later
        PredDict['scalings'].append(Scalings)
        PredDict['seqlen'].append(TStep)
        PredDict['trainfeats'].append(TrainFeatCol)
        PredDict['targfeats'].append(TargFeatCol)
        
        # Create temporal sequences of data based on daily timesteps
        X, y, TrainInd = CreateSequences(TrainFeat, TargFeat, TStep)
        
        # Separate test and train data and add to prediction dict (can't stratify when y is multicolumn)
        X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=ValidSize, random_state=0)#, stratify=y)
        PredDict['validsize'].append(ValidSize)
        PredDict['X_val'].append(X_val)
        PredDict['y_val'].append(y_val)
        
        # Use SMOTE for oversampling when dealing with imbalanced classification
        if UseSMOTE is True:
            smote = SMOTE()
            X_train_smote, y_train_smote = smote.fit_resample(X_train.reshape(X_train.shape[0], -1), y_train)
            X_train_smote = X_train_smote.reshape(X_train_smote.shape[0], X_train.shape[1], X_train.shape[2])
            # set back onto original sequenced features/labels
            PredDict['X_train'].append(X_train_smote)
            PredDict['y_train'].append(y_train_smote)
        else: # just use unsmoted sequenced data from above
            PredDict['X_train'].append(X_train)
            PredDict['y_train'].append(y_train)
            
    return PredDict, VarDFDay_scaled, VarDFDayTest_scaled


def CompileRNN(PredDict, epochNums, batchSizes, denseLayers, dropoutRt, learnRt, hiddenLscale, LossFn='mse', DynamicLR=False):
    """
    Compile the NN using the settings and data stored in the NN dictionary.
    FM Sept 2024

    Parameters
    ----------
    PredDict : dict
        Dictionary to store all the NN model metadata.
    epochNums : list of int, hyperparameter
        List of different numbers of times a dataset passes through model in training.
    batchSizes : list of int, hyperparameter
        List of different numbers of samples to work through before internal model parameters are updated.
    denseLayers : list of int, hyperparameter
        Number of neurons in fully connected layer.
    dropoutRt : list of float, hyperparameter
        Rate at which to discard data randomly to avoid overfitting.
    learnRt : list of float, hyperparameter
        Magnitude that weights adjust by to reach minimum loss.
    hiddenLscale : list of int, hyperparameter
        Number of LSTM cells based on N_samples/ k * (N_inputs + N_outputs).
    LossFn : str, optional
        Function used to calculate model loss (error), from 'mse', 'CostSensitive'
        or 'Shoreshop'. The default is 'mse'.
    DynamicLR : bool, optional
        Option for including a dynamic learning rate in the model. The default is False.

    Returns
    -------
    PredDict : dict
        Dictionary to store all the NN model metadata, now with compiled models added.

    """
    
    for mlabel in PredDict['mlabel']:
        # Index of model setup
        mID = PredDict['mlabel'].index(mlabel)
        
        # Append hyperparameters to prediction dict
        PredDict['epochN'].append(epochNums[mID])
        PredDict['batchS'].append(batchSizes[mID])
        PredDict['denselayers'].append(denseLayers[mID])
        PredDict['dropoutRt'].append(dropoutRt[mID])
        PredDict['learnRt'].append(learnRt[mID])
        PredDict['hiddenLscale'].append(hiddenLscale[mID])
        
        # inshape = (N_timesteps, N_features)
        inshape = (PredDict['X_train'][mID].shape[0], PredDict['X_train'][mID].shape[2])
        
        # GRU Model (3-layer)
        # Model = Sequential([
        #                        Input(shape=inshape), 
        #                        GRU(64, return_sequences=True),
        #                        Dropout(0.2),
        #                        GRU(64, return_sequences=True),
        #                        Dropout(0.2),
        #                        GRU(32),
        #                        Dropout(0.2),
        #                        Dense(1, activation='sigmoid')
        #                        ])
        
        # Number  of hidden layers can be decided by rule of thumb:
            # N_hidden = N_trainingsamples / (scaling * (N_input + N_output))
        N_out = len(PredDict['targfeats'][mID])
        N_hidden = round(inshape[0] / (PredDict['hiddenLscale'][mID] * (inshape[1] + N_out)))
        
        # LSTM (1 layer)
        # Input() takes input shape, used for sequential models
        # LSTM() has dimension of (batchsize, timesteps, units) and only retains final timestep (return_sequences=False)
        # Dropout() randomly sets inputs to 0 during training to prevent overfitting
        # Dense() transforms output into 2 metrics (VE and WL)
        Model = Sequential([
                            Input(shape=inshape),
                            LSTM(units=N_hidden, activation='tanh', return_sequences=False),
                            Dense(PredDict['denselayers'][mID], activation='relu'),
                            Dropout(PredDict['dropoutRt'][mID]),
                            Dense(N_out) 
                            ])
        
        # Compile model and define loss function and metrics
        if DynamicLR:
            print('Using dynamic learning rate...')
            # Use a dynamic (exponentially decreasing) learning rate
            LRschedule = ExponentialDecay(initial_learning_rate=0.001, decay_steps=10000, decay_rate=0.0001)
            Opt = tf.keras.optimizers.Adam(learning_rate=LRschedule)
        else:
            Opt = Adam(learning_rate=float(PredDict['learnRt'][mID]))
                       
        if LossFn == 'CostSensitive':
            print('Using cost sensitive loss function...')
            # Define values for false +ve and -ve and create matrix
            falsepos_cost = 1   # Inconvenience of incorrect classification
            falseneg_cost = 100 # Risk to infrastructure by incorrect classification
            binary_thresh = 0.5
            LossFn = CostSensitiveLoss(falsepos_cost, falseneg_cost, binary_thresh)
        
        elif LossFn == 'Shoreshop':
            LossFn = ShoreshopLoss
        else:
            # Just use MSE loss fn and static learning rates
            LossFn = 'mse'
        
        Model.compile(optimizer=Opt, 
                         loss=LossFn, 
                         metrics=['accuracy','mse'])
        # Save model infrastructure to dictionary of model sruns
        PredDict['model'].append(Model)
    
    return PredDict
  

def TrainRNN(PredDict, filepath, sitename, EarlyStop=False, Looped=False):
    """
    Train the compiled NN based on the training data set aside for it. Results
    are written to PredDict which is saved to a pickle file. If TensorBoard is
    used as the callback, the training history is also written to log files for
    viewing within a TensorBoard dashboard.
    FM Sept 2024

    Parameters
    ----------
    PredDict : dict
        Dictionary to store all the NN model metadata.
    filepath : str
        Filepath to save the PredDict dictionary to (for reading the trained model back in).
    sitename : str
        Name of the site of interest.
    EarlyStop : bool, optional
        Flag to include early stopping to avoid overfitting. The default is False.

    Returns
    -------
    PredDict : dict
        Dictionary to store all the NN model metadata, now with trained NN models.

    """
    if Looped is True:
        predictpath = os.path.join(filepath, sitename,'predictions','fullrun')
    else:
        predictpath = os.path.join(filepath, sitename,'predictions')
    
    if os.path.isdir(predictpath) is False:
        os.mkdir(predictpath)
    tuningpath = os.path.join(predictpath,'tuning')
    if os.path.isdir(tuningpath) is False:
        os.mkdir(tuningpath)
    filedt = dt.datetime.now().strftime('%Y%m%d-%H%M%S')
    tuningdir = os.path.join(tuningpath, filedt)
    if os.path.isdir(tuningdir) is False:
        os.mkdir(tuningdir)
            
    # Define hyperparameters to log
    HP_EPOCHS = hp.HParam('epochs', hp.Discrete(PredDict['epochN']))
    HP_BATCH_SIZE = hp.HParam('batch_size', hp.Discrete(PredDict['batchS']))
    HP_DENSE_LAYERS = hp.HParam('dense_layers', hp.Discrete(PredDict['denselayers']))
    HP_DROPOUT = hp.HParam('dropout', hp.Discrete(PredDict['dropoutRt']))
    HP_LEARNRT = hp.HParam('learn_rate', hp.Discrete(PredDict['learnRt']))

    metric_loss = 'val_loss'
    metric_accuracy = 'val_accuracy'
    
    for MLabel in PredDict['mlabel']:
        print(f"Run: {MLabel}")
        # Index of model setup
        mID = PredDict['mlabel'].index(MLabel)
        Model = PredDict['model'][mID]
    
        # TensorBoard directory to save log files to
        rundir = os.path.join(tuningdir, MLabel)
        
        HPWriter = tf.summary.create_file_writer(rundir)

        if EarlyStop:
            # Implement early stopping to avoid overfitting
            ModelCallbacks = [EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True),
                              TensorBoard(log_dir=rundir)]
        else:
            # If no early stopping, just send log data to TensorBoard for analysis
            ModelCallbacks = [TensorBoard(log_dir=rundir)]
            
        start=time.time() # start timer
        # Train the model on the training data, writing results to TensorBoard
        History = Model.fit(PredDict['X_train'][mID], PredDict['y_train'][mID], 
                            epochs=PredDict['epochN'][mID], batch_size=PredDict['batchS'][mID],
                            validation_data=(PredDict['X_val'][mID],PredDict['y_val'][mID]),
                            callbacks=[ModelCallbacks],
                            verbose=0)
        end=time.time() # end time
        
        # Write evaluation metrics to TensorBoard for hyperparam tuning
        FinalLoss, FinalAccuracy, FinalMAE = Model.evaluate(PredDict['X_val'][mID], PredDict['y_val'][mID],verbose=0)

        PredDict['loss'].append(FinalLoss)
        PredDict['accuracy'].append(FinalAccuracy)
        
        print(f"Accuracy: {FinalAccuracy}")
        # Time taken to train model
        PredDict['train_time'].append(end-start)
        print(f"Train time: {end-start} seconds")
        
        PredDict['history'].append(History)
        
        # Log hyperparameters and metrics to TensorBoard
        with HPWriter.as_default():
            hp.hparams({
                HP_EPOCHS: PredDict['epochN'][mID],
                HP_BATCH_SIZE: PredDict['batchS'][mID],
                HP_DENSE_LAYERS:PredDict['denselayers'][mID],
                HP_DROPOUT: PredDict['dropoutRt'][mID],
                HP_LEARNRT: PredDict['learnRt'][mID]
            })
            tf.summary.scalar(metric_loss, FinalLoss, step=1)
            tf.summary.scalar(metric_accuracy, FinalAccuracy, step=1)
    
    # Save trained models in dictionary for posterity
    pklpath = os.path.join(predictpath, f"{filedt+'_'+'_'.join(PredDict['mlabel'])}.pkl")
    with open(pklpath, 'wb') as f:
        pickle.dump(PredDict, f)
            
    return PredDict


def FeatImportance(PredDict, mID=0):
    """
    Calculate feature importance from trained model, using Integrated Gradients.
    FM Feb 2025
    
    Parameters
    ----------
    PredDict : dict
        Dictionary to store all the NN model metadata, now with trained NN models.
    mID : int, optional
        ID of the chosen model run stored in PredDict. The default is 0.

    Returns
    -------
    IntGradAttr : array
        Array of integrated gradient values for chosen sequence of training features.

    """
    
    def IntGrads(model, baseline, input_sample, steps=50):
        """Computes Integrated Gradients for a single input sample."""
        
        alphas = tf.linspace(0.0, 1.0, steps)  # Interpolation steps
    
        # Compute the interpolated inputs
        interp_inputs = baseline + alphas[:, tf.newaxis, tf.newaxis] * (input_sample - baseline)
        # Compute gradients at each interpolation step
        with tf.GradientTape() as tape:
            tape.watch(interp_inputs)
            preds = model(interp_inputs)
        gradients = tape.gradient(preds, interp_inputs)
    
        # Average the gradients and multiply by input difference
        avg_gradients = tf.reduce_mean(gradients, axis=0)
        return (input_sample - baseline) * avg_gradients
    
    # Define the variables
    Model = PredDict['model'][mID]
    # validation data used because goal is to explain what model learned, not robustness of model (yet)
    X_val = PredDict['X_val'][mID]

    # Define an all-zero baseline to start at
    Baseline = tf.zeros_like(X_val[-1:], dtype=tf.float32)
    # Feed the most recent validation data to test
    InSample = tf.convert_to_tensor(X_val[-1:], dtype=tf.float32)
    
    # Compute Integrated Gradients
    IntGradAttr = IntGrads(Model, Baseline, InSample)
    # Convert to NumPy for interpretation
    IntGradAttr = IntGradAttr.numpy()

    return IntGradAttr


def SHAPTest(PredDict, mID):
    """
    IN DEVELOPMENT
    FM Feb 2025

    Parameters
    ----------
    PredDict : dict
        Dictionary to store all the NN model metadata, now with trained NN models.
    mID : int
        ID of the chosen model run stored in PredDict.

    Returns
    -------
    None.

    """
    from tensorflow import keras
    
    
    model = keras.models.load_model(PredDict['model'][mID])
    X_train = PredDict['X_train'][mID]
    X_test = PredDict['X_val'][mID]

    shap.initjs()
    
    # # Select a smaller subset of training data for SHAP
    # sample_idx = np.random.choice(X_train.shape[0], 100, replace=False)  # 100 random samples
    # X_train_sample = X_train[sample_idx]  

    # Use SHAP's DeepExplainer with the function-wrapped model
    explainer = shap.GradientExplainer(model, X_train)
    
    # Compute SHAP values for the test set
    shap_values = explainer.shap_values(X_test)
    
    # Reshape X_test for SHAP summary plot (flatten time steps)
    shap.summary_plot(shap_values[0], X_test.reshape((607, -1)))


    # return shap_values


def TrainRNN_Optuna(PredDict, mlabel):
    """
    Train recurrent neural network using the package Optuna (https://optuna.readthedocs.io/en/stable/index.html). 
    Alternative to TrainRNN(). Results are written to PredDict which is saved to 
    a pickle file. If TensorBoard is used as the callback, the training history
    is also written to log files for viewing within a TensorBoard dashboard.
    FM Dec 2024

    Parameters
    ----------
    PredDict : dict
        Dictionary to store all the NN model metadata.
    mlabel : str
        Unique model name.

    Raises
    ------
    optuna
        Exception raised to tell the trainer that the trial was pruned (if hyperparams were already tested).

    Returns
    -------
    study : optuna.study.Study
        Optuna optimisation task i.e. a set of hyperparameter trials.

    """
    # Custom RMSE metric function
    def root_mean_squared_error(y_true, y_pred):
        return K.sqrt(K.mean(K.square(y_pred - y_true)))

    def set_seed(seed):
        tf.random.set_seed(seed)
        os.environ['PYTHONHASHSEED'] = str(seed)
        np.random.seed(seed)
        random.seed(seed)

    # Index of model setup
    mID = PredDict['mlabel'].index(mlabel)
    # inshape = (N_timesteps, N_features)
    inshape = (PredDict['X_train'][mID].shape[0], PredDict['X_train'][mID].shape[2])    

    # Define the LSTM model (Input No of layers should be atleast 2, including Dense layer, the inputs are initial values only)
    def CreateModel(learning_rate, batch_size, num_layers, num_nodes, dropout_rate, epochs, xtrain, ytrain, xtest, ytest):
        set_seed(42)
        min_delta = 0.001

        model = Sequential()
        model.add(Input(inshape))
        model.add(LSTM(num_nodes, activation='relu', return_sequences=True))
        model.add(Dropout(dropout_rate))
        
        for _ in range(num_layers-2):
            model.add(LSTM(num_nodes, return_sequences=True))
            model.add(Dropout(dropout_rate))
            
        model.add(LSTM(num_nodes))
        model.add(Dropout(dropout_rate)) 
        model.add(Dense(2))
                
        model.compile(loss='mse', optimizer=Adam(learning_rate=learning_rate), metrics=[root_mean_squared_error])
        
        # Train model with early stopping callback
        early_stopping_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10, min_delta=min_delta, restore_best_weights=True)
        
        history = model.fit(xtrain, ytrain, epochs=epochs, batch_size=batch_size, 
                            validation_data=(xtest, ytest), verbose=1, shuffle = False, callbacks=[early_stopping_callback])
    
        return model, history

    # Set to store unique hyperparameter configurations
    hyperparameter_set = set()

    # For Optuna, it is required to specify a Objective function which it tries to minimise (finding the global(?) minima). Here it is loss (RMSE) on validation set
    def objective(trial):
        # Define hyperparameters to optimize with Optuna
        # Set seed for os env
        os.environ['PYTHONHASHSEED'] = str(42)
        # Set random seed for Python
        random.seed(42)
        # Set random seed for NumPy
        np.random.seed(42)
        # Set random seed for TensorFlow
        tf.random.set_seed(42)
        
        ## If you want to use a parameter space and pickup hyperpameter combos randomly. Bit time consuming
        # learning_rate = trial.suggest_float('learning_rate', 1e-4, 2e-3, log=True)
        # num_layers = trial.suggest_int('num_layers', 2, 4)
        # dropout_rate = trial.suggest_float('dropout_rate', 0.1, 0.4)
        # num_nodes = trial.suggest_int('num_nodes', 100, 300)
        # batch_size = trial.suggest_int('batch_size', 24, 92)
        # epochs = trial.suggest_int('epochs', 20, 100)
    
        # If you want to specify certain values of parameters and provide as a grid of hyperparameters (optuna samplers will pickup hyperpameter combos randomly from the grid). Bit time consuming
        learning_rate = trial.suggest_categorical('learning_rate', [0.001, 0.005, 0.01])
        num_layers = trial.suggest_categorical('num_layers', [2, 3, 4])
        dropout_rate = trial.suggest_categorical('dropout_rate', [0.1, 0.2, 0.4])
        num_nodes = trial.suggest_categorical('num_nodes', [50, 100, 200, 300])
        batch_size = trial.suggest_categorical('batch_size', [24, 32, 64, 92])
        epochs = trial.suggest_categorical('epochs', [50, 80, 100, 200])
        
        # Create a tuple of the hyperparameters
        hyperparameters = (learning_rate, num_layers, dropout_rate, num_nodes, batch_size, epochs)
    
        # Check if this set of hyperparameters has already been tried
        if hyperparameters in hyperparameter_set:
            raise optuna.exceptions.TrialPruned()  # Skip this trial and suggest a new set of hyperparameters
    
        # Add the hyperparameters to the set
        hyperparameter_set.add(hyperparameters)
    
        # Make a prediction using the LSTM model for the validation set
        optimized_model, history_HPO = CreateModel(learning_rate, batch_size, num_layers, num_nodes, dropout_rate, epochs, 
                                                   PredDict['X_train'][mID], PredDict['y_train'][mID], 
                                                   PredDict['X_val'][mID], PredDict['y_val'][mID])

        ypredict = optimized_model.predict(PredDict['X_val'][mID])
        val_loss = root_mean_squared_error(PredDict['y_val'][mID], ypredict)
        
        return val_loss # Returning the Objective function value. This will be used to optimise, through creating a surrogate model for the number of trial runs
        
    # Define the pruning callback
    pruner = optuna.pruners.MedianPruner()
    
    # Create Optuna study with Bayesian optimization (TPE) and trial pruning. There is Gausian sampler (and so many other sampling options too)
    study = optuna.create_study(direction='minimize', sampler=optuna.samplers.TPESampler(), pruner=pruner)
    
    # Optimize hyperparameters
    study.optimize(objective, n_trials=50)
    
    # Retrieve best hyperparameters
    print("Best Hyperparameters:", study.best_params)
    
    # Save the study object (If required. So that it can be loaded as in the next cell without needing to re-run this whole HPO)
    # filename = f'optuna_study_SLP_H_weekly-HsTp-LSTM{n_steps}.pkl'
    # joblib.dump(study, filename)

    return study



def RunsToCSV(tuningdir,outputCSV):
    """
    Export group of hyperparameter test outputs to one CSV file.
    FM Feb 2025

    Parameters
    ----------
    tuningdir : str
        Filepath to directory holding collection of folders made by TrainRNN() 
        (holding .tfevents files).
    outputCSV : str
        Filename of output CSV with all training history stored.

    Returns
    -------
    None.

    """
    AllData = []
    
    # Iterate through all subdirectories (each subdirectory corresponds to a run)
    for rundir in os.listdir(tuningdir):
        run_path = os.path.join(tuningdir, rundir, 'validation')
        if os.path.isdir(run_path):
            # Load the TensorBoard event data
            event_acc = EventAccumulator(run_path)
            event_acc.Reload()
            
            # Extract scalar data for each metric
            for tag in event_acc.Tags()["scalars"]:
                scalar_events = event_acc.Scalars(tag)
                for event in scalar_events:
                    AllData.append({
                        "run": rundir,
                        "tag": tag,
                        "step": event.step,
                        "value": event.value,
                        "wall_time": event.wall_time
                    })
    
    # Convert all collected data into a single DataFrame
    df = pd.DataFrame(AllData)
    
    outputPath = os.path.join(tuningdir, outputCSV)
    # Extract and save to CSV
    df.to_csv(outputPath, index=False)
    print(f"Data saved to {outputPath}")
    

def FuturePredict(PredDict, ForecastDF):
    """
    Make prediction of future vegetation edge and waterline positions for transect
    of choice, using forecast dataframe as forcing and model pre-trained on past 
    observations.
    FM Nov 2024

    Parameters
    ----------
    PredDict : dict
        Dictionary to store all the NN model metadata, now with trained NN models.
    ForecastDF : DataFrame
        Per-transect dataframe of future observations, with same columns as past training data.

    Returns
    -------
    FutureOutputs : dict
        Dict storing per-model dataframes of future cross-shore waterline and veg edge predictions

    """
    # Initialise for ulti-run outputs
    FutureOutputs = {'mlabel':PredDict['mlabel'],
                     'output':[]}
    
    # For each trained model/hyperparameter set in PredDict
    for MLabel in PredDict['mlabel']:
        # Index of model setup
        mID = PredDict['mlabel'].index(MLabel)
        Model = PredDict['model'][mID]
        
        # Separate out forecast features
        # ForecastFeat = ForecastDF[['WaveHs', 'WaveDir', 'WaveTp', 'WaveAlpha', 'Runups', 'Iribarrens']]
        ForecastFeat = ForecastDF[PredDict['trainfeats'][mID]]
        
        # Sequence forecast data to shape (samples, sequencelen, variables) 
        ForecastArr, _, ForecastInd = CreateSequences(X=ForecastFeat, 
                                                      time_steps=PredDict['seqlen'][mID])
        
        # Make prediction based off forecast data and trained model
        Predictions = Model.predict(ForecastArr)
        
        # Reverse scaling to get outputs back to their original scale
        VEPredict = PredDict['scalings'][mID]['VE'].inverse_transform(Predictions[:,0].reshape(-1, 1))
        WLPredict = PredDict['scalings'][mID]['WL'].inverse_transform(Predictions[:,1].reshape(-1, 1))
        
        FutureDF = pd.DataFrame(
                   {'futureVE': VEPredict.flatten(),
                    'futureWL': WLPredict.flatten()},
                   index=ForecastInd)
        # Get rid of any potential duplicate dates from combining DFs
        FutureDF = FutureDF[~FutureDF.index.duplicated()]
        FutureOutputs['output'].append(FutureDF)
        
    return FutureOutputs


def ShorelineRMSE(FutureOutputs, TransectDFTest):
    """
    Calculate root mean square error between the CoastLearn predicted shorelines
    and their unseen test counterparts. futureVE/WL is the vegetation edge/waterline
    for the whole test period, future10dVE/WL is the vegetation edge/waterline
    for just the first 10 days predicted following the end of validation.
    FM Mar 2025

    Parameters
    ----------
    FutureOutputs : dict
        Dict storing per-model dataframes of future cross-shore waterline and veg edge predictions
    TransectDFTest : DataFrame
        Per-transect dataframe of training and target features, output from InterpVEWLWv().

    Returns
    -------
    FutureOutputs : dict
        Dict storing per-model dataframes of future cross-shore waterline and veg edge predictions,
        now with RMSE info as a dict key.

    """
    # for each model run
    for mID in range(len(FutureOutputs['mlabel'])):
        # initialise ways to store error values
        RMSElist = []
        RMSEdict = {'futureVE':None, 'future10dVE':None, 'futureWL':None, 'future10dWL':None}
        Difflist = []
        Diffdict = {'VEdiff':None, 'WLdiff':None}
        for SL in ['VE', 'WL']:
            # Define actual and predicted VE and WL
            realVals = TransectDFTest[SL]
            realVals = realVals[~realVals.index.duplicated()]
            predVals = FutureOutputs['output'][mID]['future'+SL]
            predVals = predVals[~predVals.index.duplicated()]
            # Match indexes and remove NaNs from CreateSequence moving window
            ComboDF = pd.concat([realVals, predVals], axis=1)
            ComboDF.dropna(how='any', inplace=True)
            # Calculate RMSE
            RMSEdict['future'+SL] = root_mean_squared_error(ComboDF[SL], ComboDF['future'+SL])
            RMSEdict['future10d'+SL] = root_mean_squared_error(ComboDF[SL][:10], ComboDF['future'+SL][:10])
            # Calculate distance between predicted and actual
            Diffdict[SL+'diff'] = ComboDF['future'+SL] - ComboDF[SL]
        # Add dict of VE and WL RMSEs back to per-model-run list
        RMSElist.append(RMSEdict)
        Difflist.append(pd.DataFrame(Diffdict))
    FutureOutputs['rmse'] = RMSElist
    # Add distances between actual and predicted
    FutureOutputs['XshoreDiff'] = Difflist
    return FutureOutputs
    

def SaveRMSEtoSHP(filepath, sitename, TransectInterGDFWater, CoastalDF, FutureOutputs, Subtitle=''):
    """
    Save root mean square error information (from predicted vs. test shoreline
    predictions) as a shapefile of transects.
    FM Apr 2025

    Parameters
    ----------
    filepath : str
        Local path to COASTGUARD Data folder.
    sitename : str
        Name of site of interest.
    TransectInterGDFWater : GeoDataFrame
        GeoDataFrame of cross-shore transects, intersected with waterlines.
    CoastalDF : DataFrame
        DataFrame of cross-shore transects (rows) and intersected coastal 
        timeseries/metrics (columns).
    FutureOutputs : dict
        Dict storing per-model dataframes of future cross-shore waterline and veg edge predictions,
        now with RMSE info as a dict key.
    Subtitle : str, optional
        Additional filename title to differentiate between full site model runs. The default is ''.

    Returns
    -------
    CoastalGDF : GeoDataFrame
        DataFrame of cross-shore transects (rows) and intersected coastal 
        timeseries/metrics (columns), with RMSE stats as additional columns and 
        a geometry column from the Transect GDFs.

    """
    # Strip out prediction data to be saved to shapefile
    Rows = []
    for Tr, data in FutureOutputs.items():
        if data is None:
            # Transect is missing: fill with NaNs
            Row = {
                'TransectID': Tr,
                'preddate': np.nan,
                'predVE': np.nan,
                'predWL': np.nan,
                'VE_RMSE': np.nan,
                'VE_RMSE10d': np.nan,
                'WL_RMSE': np.nan,
                'WL_RMSE10d': np.nan,
            }
        else:
            output_df = data['output'][0] if data['output'] else None
            if output_df is not None:
                preddate = output_df.index.tolist()
                FirstDt = preddate[0].strftime('%Y-%m-%d')
                LastDt = preddate[-1].strftime('%Y-%m-%d')
                predVE = output_df['futureVE'].tolist()
                predWL = output_df['futureWL'].tolist()
            else:
                preddate = predVE = predWL = np.nan
            
            rmse_dict = data['rmse'][0] if data['rmse'] else {}
            
            Row = {
                'TransectID': Tr,
                'preddate': preddate,
                'predVE': predVE,
                'predWL': predWL,
                'VE_RMSE': rmse_dict.get('futureVE', np.nan),
                'VE_RMSE10d': rmse_dict.get('future10dVE', np.nan),
                'WL_RMSE': rmse_dict.get('futureWL', np.nan),
                'WL_RMSE10d': rmse_dict.get('future10dWL', np.nan),
            }
        
        Rows.append(Row)
    
    # Create the DataFrame
    FutureOutputsDF = pd.DataFrame(Rows)
    
    # Merge
    CoastalGDF = TransectInterGDFWater.merge(FutureOutputsDF, on='TransectID')
    
    # Make a copy to edit
    CoastalSHP = CoastalGDF.copy()
    
    # Reformat fields with lists to strings
    KeyName = list(CoastalSHP.select_dtypes(include='object').columns)
    
    for Key in KeyName:
        # find a real (non-empty) value to inspect type
        realInd = next(
            (i for i, j in enumerate(CoastalSHP[Key]) if j is not None and not (isinstance(j, float) and np.isnan(j))),
            None
        )
        if realInd is None:
            continue  # skip fully empty columns
    
        val = CoastalSHP[Key].iloc[realInd]
    
        if isinstance(val, list):  # list column
            if len(val) > 0 and isinstance(val[0], (np.float64, float)):  # list of floats
                # round the floats inside lists
                CoastalSHP[Key] = CoastalSHP[Key].apply(lambda x: [round(i, 2) for i in x] if isinstance(x, list) else x)
        elif isinstance(val, (np.float64, float)):  # single float column
            # round the single floats
            CoastalSHP[Key] = CoastalSHP[Key].round(2)
        
        # finally convert all entries in this column to strings
        CoastalSHP[Key] = CoastalSHP[Key].astype(str)
    
    # Now save
    CoastalSHP.to_file(
        os.path.join(filepath, sitename, 'veglines', f'{sitename}_Transects_Intersected_Future_{FirstDt}_{LastDt}{Subtitle}.shp')
    )
    
    return CoastalGDF
    

def RMSE_Stats(FutureOutputs):
    """
    Return and print summary RMSE statistics for a full-site run of CoastLearn.
    FM May 2025

    Parameters
    ----------
    FutureOutputs : dict
        Dict storing per-model dataframes of future cross-shore waterline and veg edge predictions,
        now with RMSE info as a dict key.

    Returns
    -------
    VERMSE : list
        Per-transect list of RMSE values for predicted vegetation edge.
    WLRMSE : list
        Per-transect list of RMSE values for predicted waterline.
    VERMSE10d : list
        Per-transect list of RMSE values for predicted vegetation edge (just first 10 days).
    WLRMSE10d : list
        Per-transect list of RMSE values for predicted waterline (just first 10 days).

    """
    # Loop through RMSEs to populate lists
    VERMSE = []
    WLRMSE = []
    for Tr in FutureOutputs.keys():
        VERMSE.append(FutureOutputs[Tr]['rmse'][0]['futureVE'])
        WLRMSE.append(FutureOutputs[Tr]['rmse'][0]['futureWL'])

    # Loop through RMSEs to populate lists (first 10 days of predictions)
    VERMSE10d = []
    WLRMSE10d = []
    for Tr in FutureOutputs.keys():
        VERMSE10d.append(FutureOutputs[Tr]['rmse'][0]['future10dVE'])
        WLRMSE10d.append(FutureOutputs[Tr]['rmse'][0]['future10dWL'])
      
    # Print summary stats of transect-based RMSE values (site-wide minimum, median, maximum RMSE)
    print("minimum / median / maximum RMSE for both shorelines, full period:")
    print(f"{ np.nanmin([WLRMSE,VERMSE]) } / { np.nanmedian([WLRMSE,VERMSE]) } / { np.nanmax([WLRMSE,VERMSE]) }")
    print("minimum / median / maximum RMSE for waterline, full period:")
    print(f"{ np.nanmin(WLRMSE) } / { np.nanmedian(WLRMSE) } / { np.nanmax(WLRMSE) }")
    print("minimum / median / maximum RMSE for vegetation edge, full period:")
    print(f"{ np.nanmin(VERMSE) } / { np.nanmedian(VERMSE) } / { np.nanmax(VERMSE) }")
    
    # Print summary stats of transect-based RMSE values (site-wide minimum, median, maximum RMSE)
    # (first 10 days of predictions)
    print("minimum / median / maximum RMSE for both shorelines, first 10 days:")
    print(f"{ np.nanmin([WLRMSE10d,VERMSE10d]) } / { np.nanmedian([WLRMSE10d,VERMSE10d]) } / { np.nanmax([WLRMSE10d,VERMSE10d]) }")
    print("minimum / median / maximum RMSE for waterline, first 10 days:")
    print(f"{ np.nanmin(WLRMSE10d) } / { np.nanmedian(WLRMSE10d) } / { np.nanmax(WLRMSE10d) }")
    print("minimum / median / maximum RMSE for vegetation edge, first 10 days:")
    print(f"{ np.nanmin(VERMSE10d) } / { np.nanmedian(VERMSE10d) } / { np.nanmax(VERMSE10d) }")
    
    return VERMSE, WLRMSE, VERMSE10d, WLRMSE10d


#%% CLASSIFICATION OF IMPACTS ###

def ClassifyImpact(TransectDF, Method='pcnt'):
    """
    Classify impact of coastal shoreline changes predicted by CoastLearn.
    FM Mar 2025

    Parameters
    ----------
    TransectDF : DataFrame
        Dataframe of per-transect coastal metrics/variables in timeseries.
        Usually 'FutureOutputsClean[Tr]['output'][0]' is passed.
    Method : str, optional
        Method to use for calculating the impact factor. The default is 'pcnt'. Choose from:
        'pcnt' (use the 5th and 25th percentiles of cross-shore distance 
                to classify the bounds of high and medium impact); 
        'slope' (use the 5th and 25th percentiles of cross-shore timeseries slope 
                to classify the bounds of high and medium impact); 
        'combi' (use a combination of distance and slope percentiles). 

    Returns
    -------
    ImpactClass : dict
        Dict of impact classifications for futureVE and futureWL, plus the bounds
        calculated to apply the classifications.

    """
    
    # For each type of shoreline
    if 'future' in TransectDF.columns.any(): # if classifying future shorelines
        SLkeys = ['futureWL','futureVE']
        # Define dictionary
        ImpactClass = {'futureWL':[], 'futureVE':[]}
    else:
        SLkeys = ['WL','VE']
        # Define dictionary
        ImpactClass = {'WL':[], 'VE':[]}
    for SL in SLkeys:
        # If using percentiles
        if Method=='pcnt':
            # Mid impact is between 5th and 25th percentile
            MidQuant = TransectDF[SL].quantile(q=0.25)
            # High impact is lower than 5th percentile
            LowQuant =  TransectDF[SL].quantile(q=0.05)
            ImpactClass[SL+'_MidLim'] = MidQuant
            ImpactClass[SL+'_HighLim'] = LowQuant
            
            for i in range(len(TransectDF[SL])):
                if TransectDF[SL].iloc[i] <= LowQuant:
                    # High impact
                    ImpactClass[SL].append(3)
                elif LowQuant < TransectDF[SL].iloc[i] < MidQuant:
                    # Mid impact
                    ImpactClass[SL].append(2)
                else:
                    # Low impact
                    ImpactClass[SL].append(1)
        
        elif Method=='combi':
            # Mid impact is between 5th and 25th percentile
            MidQuant = TransectDF[SL].quantile(q=0.25)
            # High impact is lower than 5th percentile
            LowQuant =  TransectDF[SL].quantile(q=0.05)
            ImpactClass[SL+'_MidLim'] = MidQuant
            ImpactClass[SL+'_HighLim'] = LowQuant
            
            Slopes = TransectDF[SL].diff() / TransectDF.index.to_series().diff().dt.total_seconds()
            # Mid impact is between 5th and 25th percentile
            MidSlope = Slopes.quantile(q=0.25)
            # High impact is lower than 5th percentile
            SteepSlope =  Slopes.quantile(q=0.05)
            ImpactClass[SL+'_MidSlope'] = MidSlope
            ImpactClass[SL+'_HighSlope'] = SteepSlope
            
            for i in range(len(TransectDF[SL])):
                # Use quantile as base class
                if TransectDF[SL].iloc[i] <= LowQuant:
                    # High impact
                    ImpactClass[SL].append(3)
                elif TransectDF[SL].iloc[i] <= MidQuant:
                    # Mid impact
                    ImpactClass[SL].append(2)
                else:
                    # Low impact
                    ImpactClass[SL].append(1)
                # Add mid and steep slopes on top
                if Slopes.iloc[i] <= SteepSlope:
                    ImpactClass[SL][-1] = max(ImpactClass[SL][-1],3)
                elif Slopes.iloc[i] <= MidSlope:
                    ImpactClass[SL][-1] = max(ImpactClass[SL][-1],2)
        
        # Or if using slope        
        else:
            Slopes = TransectDF[SL].diff() / TransectDF.index.to_series().diff().dt.total_seconds()
            # Mid impact is between 5th and 25th percentile
            MidQuant = Slopes.quantile(q=0.25)
            # High impact is lower than 5th percentile
            LowQuant =  Slopes.quantile(q=0.05)
            ImpactClass[SL+'_MidSlope'] = MidQuant
            ImpactClass[SL+'_HighSlope'] = LowQuant
            
            for i in range(len(TransectDF[SL])):
                if Slopes.iloc[i] <= LowQuant:
                    ImpactClass[SL].append(3)
                elif LowQuant < Slopes.iloc[i] < MidQuant:
                    ImpactClass[SL].append(2)
                else:
                    ImpactClass[SL].append(1)
            
        ImpactClass[SL] = pd.Series(ImpactClass[SL], index=TransectDF.index)
        
    if 'future' in TransectDF.columns.any():
        ImpactClass['Sum'] = ImpactClass['futureWL'] + ImpactClass['futureVE']
    else:
        ImpactClass['Sum'] = ImpactClass['WL'] + ImpactClass['VE']
    


    return ImpactClass


# ----------------------------------------------------------------------------------------
#%% CLUSTERING FUNCTIONS ###

def Cluster(TransectDF, ValPlots=False):
    """
    Classify coastal change indicator data into low, medium or high impact from hazards,
    using a SpectralCluster clustering routine.
    FM Sept 2024

    Parameters
    ----------
    TransectDF : DataFrame
        Dataframe of single cross-shore transect, with timeseries of satellite-derived metrics attached.
    ValPlots : bool, optional
        Plot validation plots of silhouette score and inertia. The default is False.

    Returns
    -------
    VarDFClust : DataFrame
        Dataframe of just coastal metrics/variables in timeseries, with cluster values attached to each timestep.

    """
    
    from Toolshed import PredictionsPlotting
    
    # Fill nans factoring in timesteps for interpolation
    TransectDF.replace([np.inf, -np.inf], np.nan, inplace=True)
    VarDF = TransectDF.interpolate(method='time', axis=0)
    
    # VarDF = VarDF[['distances', 'wlcorrdist','TZwidth','WaveHs']]
    
    VarDF_scaled = StandardScaler().fit_transform(VarDF)
    
    # Apply PCA to reduce the dimensions to 3D for visualization
    pca = PCA(n_components=2)
    pca_VarDF = pca.fit_transform(VarDF_scaled)
    eigenvectors = pca.components_
    variances = pca.explained_variance_ratio_
    
    # ClusterMods = {'':SpectralClustering(n_clusters=3, eigen_solver='arpack', random_state=42)}
    # for Mod in ClusterMods.keys():
        
    ClusterMod = SpectralClustering(n_clusters=3, 
                                    eigen_solver='arpack',
                                    n_components=len(VarDF.columns), 
                                    random_state=42)
    # Map labels to cluster IDs based on cluster centres and their distance to eigenvectors
    ClusterMod.fit(VarDF_scaled)
    VarDF['Cluster'] = ClusterMod.labels_
    # 
    # ClusterCentres = np.array([pca_VarDF[VarDF['Cluster'] == i].mean(axis=0) for i in range(3)])
    
    # # Define cluster labels using identified centres
    # HighImpact = np.argmax(ClusterCentres[:, 0])
    # LowImpact = np.argmax(ClusterCentres[:, 1])
    # MediumImpact = (set([0,1,2]) - {HighImpact, LowImpact}).pop()
    # # Map labels to cluster IDs
    # ClusterToImpact = {HighImpact:'High',
    #                    MediumImpact:'Medium',
    #                    LowImpact:'Low'}
    # ImpactLabels = [ClusterToImpact[Cluster] for Cluster in VarDF['Cluster']]
    # VarDFClust = VarDF.copy()
    # VarDFClust['Impact'] = ImpactLabels

    # Create a DataFrame for PCA results and add cluster labels
    pca_df = pd.DataFrame(data=pca_VarDF, columns=['PC1', 'PC2'])
    pca_df['Cluster'] = ClusterMod.labels_

    # Visualization of clusters
    # Example clustered timeseries using one or two variables
    if ValPlots is True:
        PredictionsPlotting.PlotClusteredTS(VarDF)

    scale_factor = 0.5
    # Plot the clusters in the PCA space
    PredictionsPlotting.PlotCluster(VarDF, pca_df, scale_factor, eigenvectors, variances)
        
        
    return VarDF


def ClusterKMeans(TransectDF, ValPlots=False):
    """
    Classify coastal change indicator data into low, medium or high impact from hazards,
    using a KMeans clustering routine.
    FM Sept 2024

    Parameters
    ----------
    TransectDF : DataFrame
        Dataframe of single cross-shore transect, with timeseries of satellite-derived metrics attached.
    ValPlots : bool, optional
        Plot validation plots of silhouette score and inertia. The default is False.

    Returns
    -------
    VarDFClust : DataFrame
        Dataframe of just coastal metrics/variables in timeseries, with cluster values attached to each timestep.

    """
    from Toolshed import PredictionsPlotting
    
    # Define variables dataframe from transect dataframe by removing dates and transposing
    VarDF = TransectDF.drop(columns=['TransectID', 'dates'])
    VarDF.interpolate(method='nearest', axis=0, inplace=True) # fill nans using nearest
    VarDF.interpolate(method='linear', axis=0, inplace=True) # if any nans left over at start or end, fill with linear
    VarDF_scaled = StandardScaler().fit_transform(VarDF)
    
    
    # Fit k-means clustering to data iteratively over different cluster sizes
    k_n = range(2,15)
    # Inertia = compactness of clusters i.e. total variance within a cluster
    # Silhouette score = how similar object is to its own cluster vs other clusters 
    inertia = []
    sil_scores = []
    
    for k in k_n:
        kmeansmod = KMeans(n_clusters=k, random_state=42)
        kmeansmod.fit(VarDF_scaled)
        inertia.append(kmeansmod.inertia_)
        sil_scores.append(silhouette_score(VarDF_scaled, kmeansmod.labels_))
    
    # Apply PCA to reduce the dimensions to 3D for visualization
    pca = PCA(n_components=3)
    pca_VarDF = pca.fit_transform(VarDF_scaled)
    eigenvectors = pca.components_

    # Create a DataFrame for PCA results and add cluster labels
    pca_df = pd.DataFrame(data=pca_VarDF, columns=['PC1', 'PC2', 'PC3'])
    pca_df['Cluster'] = kmeansmod.labels_
    
    
    if ValPlots is True:
        PredictionsPlotting.PlotClusterElbowSil(k_n, inertia, sil_scores)
    
    
    # # Fit the KMeans model with the chosen number of clusters
    # # Clusters are informed by 'impact' levels low, medium and high
    # optimal_k = 3
    # tic = timeit.default_timer() # start timer
    # kmeansmod = KMeans(n_clusters=optimal_k, random_state=42)
    # kmeansmod.fit(VarDF_scaled)
    # toc = timeit.default_timer() # stop timer
    
    # # Analyze the clustering results
    # VarDF['Cluster'] = kmeansmod.labels_
    
    ClusterMods = {'kmeans':KMeans(n_clusters=3, random_state=42),
                   'spectral':SpectralClustering(n_clusters=3, eigen_solver='arpack', random_state=42)}
    for Mod in ClusterMods.keys():
        
        ClusterMods[Mod].fit(VarDF_scaled)
        VarDF[Mod+'Cluster'] = ClusterMods[Mod].labels_
        ClusterMeans = VarDF.groupby(Mod+'Cluster').mean()
        
        ClusterCentres = np.array([pca_VarDF[VarDF[Mod+'Cluster'] == i].mean(axis=0) for i in range(3)])

        HighImpact = np.argmax(ClusterCentres[:, 0])
        LowImpact = np.argmax(ClusterCentres[:, 1])
        MediumImpact = (set([0,1,2]) - {HighImpact, LowImpact}).pop()
        
        ClusterToImpact = {HighImpact:'High',
                           MediumImpact:'Medium',
                           LowImpact:'Low'}
        ImpactLabels = [ClusterToImpact[Cluster] for Cluster in VarDF[Mod+'Cluster']]
        VarDFClust = VarDF.copy()
        VarDFClust[Mod+'Impact'] = ImpactLabels
        
        # HighImpact = ClusterMeans[(ClusterMeans['distances'] == ClusterMeans['distances'].min()) & # landward VE
        #                           (ClusterMeans['wlcorrdist'] == ClusterMeans['wlcorrdist'].min()) & # landward WL
        #                           (ClusterMeans['waterelev'] == ClusterMeans['waterelev'].max()) & # high water
        #                           (ClusterMeans['beachwidth'] == ClusterMeans['beachwidth'].min()) & # narrow width
        #                           (ClusterMeans['TZwidth'] == ClusterMeans['TZwidth'].min()) & # narrow TZ
        #                           (ClusterMeans['WaveHs'] == ClusterMeans['WaveHs'].max()) & # high waves
        #                           (ClusterMeans['WaveTp'] == ClusterMeans['WaveTp'].max())].index[0] # long period
        
        # LowImpact = ClusterMeans[(ClusterMeans['distances'] == ClusterMeans['distances'].max()) & # seaward VE
        #                           (ClusterMeans['wlcorrdist'] == ClusterMeans['wlcorrdist'].max()) & # seaward WL
        #                           (ClusterMeans['waterelev'] == ClusterMeans['waterelev'].min()) & # low water
        #                           (ClusterMeans['beachwidth'] == ClusterMeans['beachwidth'].max()) & # wide width
        #                           (ClusterMeans['TZwidth'] == ClusterMeans['TZwidth'].max()) & # wide TZ
        #                           (ClusterMeans['WaveHs'] == ClusterMeans['WaveHs'].min()) & # low waves
        #                           (ClusterMeans['WaveTp'] == ClusterMeans['WaveTp'].min())].index[0] # short period
        # AllClusters = set([0,1,2])
        # MediumImpact = (AllClusters - set([HighImpact, LowImpact])).pop()

        # Cluster to impact
        # ClusterToImpact = {'High': HighImpact,
        #                    'Medium':MediumImpact,
        #                    'Low':LowImpact}
        # VarDF['Impact'] = VarDF[Mod+'Cluster'].map(ClusterToImpact)
        
        # inertia.append(ClusterMods[Mod].inertia_)
        # sil_scores.append(silhouette_score(VarDF_scaled, ClusterMods[Mod].labels_))
    
        # Create a DataFrame for PCA results and add cluster labels
        pca_df = pd.DataFrame(data=pca_VarDF, columns=['PC1', 'PC2', 'PC3'])
        pca_df['Cluster'] = ClusterMods[Mod].labels_
    
        # Optional: Visualization of clusters
        if ValPlots is True:
            PredictionsPlotting.PlotClusterVisuals(VarDF, Mod, pca_df, eigenvectors)
        
        
    return VarDFClust
