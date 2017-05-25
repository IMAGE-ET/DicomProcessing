import SimpleITK as sitk
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import skimage.io as io
import sys,os,fnmatch,shutil
import collections,random
from ipywidgets import interact, interactive, fixed, interact_manual
import dicom
csvFILE='./deepmedic/examples/CTdataset/coords_lm.csv'
serverPATH='./deepmedic/examples/CTdataset/Dataset/'
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


def searchPath(reader,input_dir):
    series=[]
    paths=[]
    tags=[]
    nbdl=1


    for cur_dir, subdirs, files in os.walk( input_dir ):
        series_found = reader.GetGDCMSeriesIDs( cur_dir )
        if series_found:
            series.append( series_found )
            paths.append( cur_dir )
            tag = "_".join( cur_dir.split( os.path.sep )[ nbdl-1 : ] ).replace(" ", "_")
            tags.append( tag )
    return series,paths,tags

################### Parse CSV file with 3D landmarks ####################
land3D,landmarkList=parseLandmarks(csvFILE)
print('Reading 3D landmark coordinates',csvFILE,'...')

reader = sitk.ImageSeriesReader()
series,paths,tags=searchPath(reader,serverPATH)

#Process series: segment face, extract contour, compute mesh
for ii, path in enumerate( paths ):
    for jj, serie in enumerate(series[ii]):
        print (path,serie)
        try:
            dicom_names = reader.GetGDCMSeriesFileNames( path, serie)
        except:
            continue
        dcmf = dicom.read_file(dicom_names[0])
        serie=dcmf.data_element("SeriesInstanceUID").value
        seriename=dcmf.data_element("SeriesDescription").value
        seriename=seriename.replace(" ","").replace(".","")
        seriesnumber=dcmf.data_element("SeriesNumber").value
        pid=dcmf.data_element("PatientID").value

        patientdir=pid
        seriedir=str(seriename)+"_"+str(seriesnumber)

        if(land3D[patientdir]):
            print(serie,patientdir,seriedir)
            #saving image in a predefined path style
            # newdir="./"+str(patientdir)+"/"+seriedir+"/"
            # newfiledir=newdir+"newname.dcm"
            # if not os.path.exists(newfiledir):
            #     os.makedirs(newdir)
            #     dcmf.save_as(newfiledir)
            # else:
            #     print("FILE exists already...check the script","...")
            #     programPause = raw_input("Press the <ENTER> key to continue...")
