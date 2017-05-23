import SimpleITK as sitk
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import skimage.io as io
import sys,os,fnmatch,shutil
import collections,random

# Collect all DICOM series
reader = sitk.ImageSeriesReader()
for cur_dir, subdirs, files in os.walk( base_directory ):
    series_found = reader.GetGDCMSeriesIDs( cur_dir )
    if series_found:
        series.append( series_found )
        paths.append( cur_dir )
        tag = "_".join( cur_dir.split( os.path.sep )[ nbdl-1 : ] ).replace(" ", "_")
        tags.append( tag )

#Process series: segment face, extract contour, compute mesh
for ii, path in enumerate( paths ):
    for jj, serie in enumerate(series[ii]):
        print path
        try:
            dicom_names = reader.GetGDCMSeriesFileNames( path, serie)
        except:
            continue
        if len(dicom_names):
            reader.SetFileNames(dicom_names)
            try:
                image = reader.Execute()
            except:
                continue
                
                
                
# data_directory = os.path.dirname(fdata("CIRS057A_MR_CT_DICOM/readme.txt"))
            # Global variable 'selected_series' is updated by the interact function
# selected_series = ''
# def DICOM_series_dropdown_callback(series_to_load, series_dictionary):
    # global selected_series
               # Print some information about the series from the meta-data dictionary
              # DICOM standard part 6, Data Dictionary: http://medical.nema.org/medical/dicom/current/output/pdf/part06.pdf
    # img = sitk.ReadImage(series_dictionary[series_to_load][0])
    # tags_to_print = {'0010|0010': 'Patient name: ', 
                     # '0008|0060' : 'Modality: ',
                     # '0008|0021' : 'Series date: ',
                     # '0008|0080' : 'Institution name: ',
                     # '0008|1050' : 'Performing physician\'s name: '}
    # for tag in tags_to_print:
        # try:
            # print(tags_to_print[tag] + img.GetMetaData(tag))
        # except: # Ignore if the tag isn't in the dictionary
            # pass
    # selected_series = series_to_load                    

             # Directory contains multiple DICOM studies/series, store
             # in dictionary with key being the seriesID
# reader = sitk.ImageSeriesReader()
# series_file_names = {}
# series_IDs = reader.GetGDCMSeriesIDs(data_directory)
            # Check that we have at least one series
# if series_IDs:
    # for series in series_IDs:
        # series_file_names[series] = reader.GetGDCMSeriesFileNames(data_directory, series)
    
    # interact(DICOM_series_dropdown_callback, series_to_load=series_IDs, series_dictionary=fixed(series_file_names)); 
# else:
    # print('Data directory does not contain any DICOM series.')