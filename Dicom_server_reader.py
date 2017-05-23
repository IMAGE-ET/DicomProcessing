import SimpleITK as sitk
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import skimage.io as io
import sys,os,fnmatch,shutil
import collections,random

csvFILE='./coords_lm.csv'

#################### 3D landmark parsing method ####################
def parseLandmarks(filename):
    file=open(filename)
    count=0
    landmarkDict= collections.defaultdict(dict)
    landmarkList= []
    for line in file:
        count=count+1
        line=line.replace('"','')
        tokens=line.strip("\n").split(",")
        name=tokens[0].strip("\t ^")
        if count==1:
            #Parse Landmark names
            for i in range(7,64,3):
                tk=tokens[i].split(".")
                landmark=tk[0]
                coord=tk[1]
                if(len(tk)==3):
                    ori=tk[2]
                    landmark=landmark+'.'+ori
                else:
                    ori=''
                    landmark=landmark
                landmarkList.append(landmark)
        else:
            listindex=0
            for i in range(7,64,3):
                lx=tokens[i]
                ly=tokens[i+1]
                lz=tokens[i+2]
                #if(lx!='' and ly!='' and lz!=''):
                #    print(lx,ly,lz)
                #use the previous if in case we dont want to save blank landmarks in the dictionary
                try:
                    landmarkDict[name][landmarkList[listindex]]=[lx,ly,lz]
                except:
                    print(name,'ERROR adding name')
                listindex=listindex+1
    print('Number of 3D Landmarks',len(landmarkList))
    print('Number of Subjects',len(landmarkDict))
    return landmarkDict,landmarkList
    file.close()
    
####################  Path reading methods for DCM images:  ####################
def searchPath(input_dir):
    matches = []
    dirs = []
    for root, dirnames, filenames in os.walk(input_dir):
        for filename in fnmatch.filter(filenames, '*.dcm'):
            matches.append(os.path.join(root, filename));
            if root not in dirs:
                dirs.append(root)
    return (matches, dirs)

def readDCMdir(dirs):
    isr = sitk.ImageSeriesReader()
    seriessets = []
    for d in dirs:
        series = isr.GetGDCMSeriesIDs(d)
        for s in series:
            files = isr.GetGDCMSeriesFileNames(d, s)
            print s, d, len(files)
            seriessets.append([s, d, files])
    return seriessets
    
    
################### Parse CSV file with 3D landmarks ####################
land3D,landmarkList=parseLandmarks(csvFILE)
print('Reading 3D landmark coordinates',csvFILE,'...')



    
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
    
    
    # def writeSeries(imagein,outDir):
    # writer = sitk.ImageSeriesWriter()
    # filenames = [ outDir+'MR_{0:04}.dcm'.format(i) for i in
# range(imagein.GetSize()[2])]
    # writer.SetFileNames(filenames)
    # writer.Execute(imagein)
    
    
    
 ##### If using pydicom

# import dicom

# dcmf = dicom.read_file("rtdcmf.dcm")
# print(dcmf)
# dcmf.PatientName
# name=dcmf.data_element("PatientsName")
# print(name)
# dcmf.PatientID
# id==dcmf.data_element("PatientsID")
# dcmf.SeriesNumber
# dcmf.dir("pat")
# tag= dcmf.PatientSetupSequence[0]
# print(tag)
# dcmf.save_as("newname.dcm")
            
    
    
    
    