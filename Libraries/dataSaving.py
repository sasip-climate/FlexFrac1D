import numpy as np
import pandas as pd

from os import path

projectDirectory = "C:/Users/tlilia/Documents/ProjetSASIP/"
databaseDirectory = 'database/'
defaultParameters = {'h':[], 'wvLength':[], 'wvAmpl':[], 'initialLength':[], 'E':[], 'nu':[], 'G':[], 'multiFrac':[], 'EnergyType':[], 'DispType':[]}

def createAndLoadDataFrame(dataFileName, parameters, results=['xc', 'FloeSizes']):
    ''' Returns DataFrame from cross product of parameter values and from already computed data in pickle dataFileName
    Inputs:
        dataFileName (string): name of the pickle file containing existing data
        parameters (dic): dictionary of parameters name and their values as a list or np.array
        results (list of str): names of computed 'outputs' 
    Outputs:
        df (pd.DataFrame): data consisting of pm1|pm2|...|pm?||out1|out2
    '''
    assert ('FloeSizes' in results) and ('xc' in results)
    assert parameters.keys() == defaultParameters.keys()
    
    
    df = None
    # Creates columns for the parameters (cross product)
    for key, value in parameters.items():
        if df is None:
            df = pd.DataFrame(data=value, columns=[key])
        else:
            parameterValues = pd.Series(data=value, name=key)
            df = df.merge(parameterValues, how='cross', suffixes=('', '_2'))
    
    # Creates columns for results (xc and FloeSizes for now)
    outputs = {res: [np.NaN]*df.shape[0] for res in results}
    outputs['computed'] = [False]*df.shape[0]
    df = df.assign(**outputs)

    # check the data file for already computed values
    if path.isfile(dataFileName) and path.getsize(dataFileName) > 0:
        oldDf = pd.read_pickle(dataFileName)

        # if same columns, merge on parameters to not duplicate already computed values
        if np.array_equal(df.columns, oldDf.columns):
            # replace the data from oldDf
            df = oldDf.merge(df, how='outer', on=list(parameters.keys()), suffixes=('', '_2'))
            
            # update if computed
            df['computed'] = df['computed'].fillna(False) | df['computed_2'].fillna(False)
            
            # Drop the _2 since every output is NaN at this moment
            for col in results+['computed']:
                if col+'_2' in df.columns:
                    df.drop(columns=col+'_2', inplace=True) 

            
        # else is harder -> for now we just overwrite it...
        # TODO: we want to be able to add a new parameter or a new output and insert NaN where not given/computed      
        else:
            raise ValueError('TODO in code: Unsual parameter/output given')

    return df

def loadDataFrame(dataFileName):
    ''' Loads DataFrame from pickle file dataFileName
    '''
    # check the data file for already computed values
    if path.isfile(dataFileName) and path.getsize(dataFileName) > 0:
        return pd.read_pickle(dataFileName)
    else:
        return None
    
def saveDataFrame(dataFrame, dataFileName):
    ''' Saves dataFrame in pickle file dataFileName
    '''
    dataFrame.to_pickle(dataFileName)
    
class listForPandas:
    ''' Encapsulation class to save lists and arrays in pandas.DataFrame cells without pandas unpacking them
    '''
    def __init__(self, list):
        self.content = list
    def __repr__(self):
        return str(self.content)