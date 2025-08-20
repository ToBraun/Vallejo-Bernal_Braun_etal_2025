# Copyright (C) 2025 by
# Tobias Braun <tobraun@pik-potsdam.de> 
# All rights reserved.
# GNU General Public License v3.0.

'''
Utility functions for the atmospheric river detection tool underlying the 
PIKART version 1 catalog.
'''

# %% IMPORT MODULES

# Standard library imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %% FIND SEQUENCES OF ONES IN A BINARY TIME SERIES

def rle(array):
    """
    Run-length encoding for an array.
    
        Parameters
    ----------
    array : array-like
        Input array to be encoded.
    
    Returns
    -------
    values : ndarray
        Array containing unique values in the input array.
    counts : ndarray
        Array containing run lengths for each unique value.
    """
    array = np.asarray(array)
    changes = np.where(array[:-1] != array[1:])[0] + 1
    
    indices = np.concatenate(([0], changes, [len(array)]))
    counts = np.diff(indices)
    values = array[indices[:-1]]
    
    return values, counts

# %% NUMBER OF AR EVENTS AND DURATION

def count_ARevents(mask):
    """
    Splits the AR mask into AR events and computes event counts and durations.
    
    Parameters
    ----------
    mask : ndarray
        3D array representing the AR mask.
    
    Returns
    ------- 
    a_count : ndarray
        Array containing counts of AR events at each grid point.
    a_dur : ndarray
        Array containing mean durations of AR events at each grid point.
    """
    _, nlat, nlon = mask.shape
    a_count, a_dur = np.zeros((nlat,nlon)), np.zeros((nlat,nlon))
    for i in range(nlat):
        for j in range(nlon):
            # grab mask time series
            tmp_mask = mask[:, i, j]
            # run length time encoding and get number of AR events
            tmp_bool, tmp_freqs = rle(tmp_mask)#rle(tmp_mask).T
            tmp_eventidx = np.where(tmp_bool == 1)[0]      
            # number of AR events and their mean duration
            L = tmp_eventidx.size
            DUR = np.mean(tmp_freqs[tmp_eventidx]*6)
            a_count[i,j] = L
            a_dur[i,j] = DUR
    return a_count, a_dur




# %% READ CSV RECORD

def read_csv_record(file_name):
    """
    Read a CSV file containing geographic contour and axis data and convert array-like 
    string columns to NumPy arrays.

    Parameters
    ----------
    file_name : str
        Path to the CSV file to read.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the CSV data with the following conversions:
        - Columns `contour_lon`, `contour_lat`, `axis_lon`, `axis_lat` are converted 
          from string representations of arrays (e.g., "[1 2 3]") to `numpy.ndarray` 
          of floats.
        - Column `id` is cast to `int`.
        - Column `time` is kept as string (`str`).
        - Column `deleted` is kept as string (`str`).
        - Column `is_relaxed` is cast to `bool`.
        - Column `why_relaxed` is kept as string (`str`).

    Notes
    -----
    The function expects array-like columns to be space-separated and enclosed in 
    square brackets (e.g., "[1.0 2.0 3.0]"). Other columns will be read according to 
    the specified `dtype`.
    """
    def conv_array(text):
        '''Convert array texts to ndarray'''
        text = text.replace('[','').replace(']','')
        array = np.array(text.split()).astype('float')
        return array

    conv_keys = ['contour_lon', 'contour_lat', 'axis_lon', 'axis_lat']

    converters = dict([(key_i, conv_array) for key_i in conv_keys])

    dtypes={'id': 'int', 
            'time': 'str',
            'deleted': 'str',
            'is_relaxed': 'bool',
            'why_relaxed': 'str'}

    ardf = pd.read_csv(file_name, dtype=dtypes, converters=converters)

    return ardf


    

# %% SINGULAR SPECTRUM ANALYSIS

class SSA(object):
    
    __supported_types = (pd.Series, np.ndarray, list)
    
    def __init__(self, tseries, L, save_mem=True):
        """
        Decomposes the given time series with a singular-spectrum analysis. Assumes the values of the time series are
        recorded at equal intervals.
        
        Parameters
        ----------
        tseries : The original time series, in the form of a Pandas Series, NumPy array or list. 
        L : The window length. Must be an integer 2 <= L <= N/2, where N is the length of the time series.
        save_mem : Conserve memory by not retaining the elementary matrices. Recommended for long time series with
            thousands of values. Defaults to True.
        
        Note: Even if an NumPy array or list is used for the initial time series, all time series returned will be
        in the form of a Pandas Series or DataFrame object.
        """
        
        # Tedious type-checking for the initial time series
        if not isinstance(tseries, self.__supported_types):
            raise TypeError("Unsupported time series object. Try Pandas Series, NumPy array or list.")
        
        # Checks to save us from ourselves
        self.N = len(tseries)
        if not 2 <= L <= self.N/2:
            raise ValueError("The window length must be in the interval [2, N/2].")
        
        self.L = L
        self.orig_TS = pd.Series(tseries)
        self.K = self.N - self.L + 1
        
        # Embed the time series in a trajectory matrix
        self.X = np.array([self.orig_TS.values[i:L+i] for i in range(0, self.K)]).T
        
        # Decompose the trajectory matrix
        self.U, self.Sigma, VT = np.linalg.svd(self.X)
        self.d = np.linalg.matrix_rank(self.X)
        
        self.TS_comps = np.zeros((self.N, self.d))
        
        if not save_mem:
            # Construct and save all the elementary matrices
            self.X_elem = np.array([ self.Sigma[i]*np.outer(self.U[:,i], VT[i,:]) for i in range(self.d) ])

            # Diagonally average the elementary matrices, store them as columns in array.           
            for i in range(self.d):
                X_rev = self.X_elem[i, ::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]
            
            self.V = VT.T
        else:
            # Reconstruct the elementary matrices without storing them
            for i in range(self.d):
                X_elem = self.Sigma[i]*np.outer(self.U[:,i], VT[i,:])
                X_rev = X_elem[::-1]
                self.TS_comps[:,i] = [X_rev.diagonal(j).mean() for j in range(-X_rev.shape[0]+1, X_rev.shape[1])]
            
            self.X_elem = "Re-run with save_mem=False to retain the elementary matrices."
            
            # The V array may also be very large under these circumstances, so we won't keep it.
            self.V = "Re-run with save_mem=False to retain the V matrix."
        
        # Calculate the w-correlation matrix.
        self.calc_wcorr()
        
    def return_eigenvalues(self):
        return self.Sigma**2
            
    def components_to_df(self, n=0):
        """
        Returns all the time series components in a single Pandas DataFrame object.
        """
        if n > 0:
            n = min(n, self.d)
        else:
            n = self.d
        
        # Create list of columns - call them F0, F1, F2, ...
        cols = ["F{}".format(i) for i in range(n)]
        return pd.DataFrame(self.TS_comps[:, :n], columns=cols, index=self.orig_TS.index)
            
    
    def reconstruct(self, indices):
        """
        Reconstructs the time series from its elementary components, using the given indices. Returns a Pandas Series
        object with the reconstructed time series.
        
        Parameters
        ----------
        indices: An integer, list of integers or slice(n,m) object, representing the elementary components to sum.
        """
        if isinstance(indices, int): indices = [indices]
        
        ts_vals = self.TS_comps[:,indices].sum(axis=1)
        return pd.Series(ts_vals, index=self.orig_TS.index)
    
    def calc_wcorr(self):
        """
        Calculates the w-correlation matrix for the time series.
        """
             
        # Calculate the weights
        w = np.array(list(np.arange(self.L)+1) + [self.L]*(self.K-self.L-1) + list(np.arange(self.L)+1)[::-1])
        
        def w_inner(F_i, F_j):
            return w.dot(F_i*F_j)
        
        # Calculated weighted norms, ||F_i||_w, then invert.
        F_wnorms = np.array([w_inner(self.TS_comps[:,i], self.TS_comps[:,i]) for i in range(self.d)])
        F_wnorms = F_wnorms**-0.5
        
        # Calculate Wcorr.
        self.Wcorr = np.identity(self.d)
        for i in range(self.d):
            for j in range(i+1,self.d):
                self.Wcorr[i,j] = abs(w_inner(self.TS_comps[:,i], self.TS_comps[:,j]) * F_wnorms[i] * F_wnorms[j])
                self.Wcorr[j,i] = self.Wcorr[i,j]
    
    def plot_wcorr(self, min=None, max=None):
        """
        Plots the w-correlation matrix for the decomposed time series.
        """
        if min is None:
            min = 0
        if max is None:
            max = self.d
        
        if self.Wcorr is None:
            self.calc_wcorr()
        
        ax = plt.imshow(self.Wcorr)
        plt.xlabel(r"$\tilde{F}_i$")
        plt.ylabel(r"$\tilde{F}_j$")
        plt.colorbar(ax.colorbar, fraction=0.045)
        ax.colorbar.set_label("$W_{i,j}$")
        plt.clim(0,1)
        
        # For plotting purposes:
        if max == self.d:
            max_rnge = self.d-1
        else:
            max_rnge = max
        
        plt.xlim(min-0.5, max_rnge+0.5)
        plt.ylim(max_rnge+0.5, min-0.5)
        