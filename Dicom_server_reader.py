import SimpleITK as sitk
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import skimage.io as io
import glob,sys,os,fnmatch,shutil,string,unicodedata
import collections,random
from ipywidgets import interact, interactive, fixed, interact_manual
import dicom
from sets import Set

csvFILE='./deepmedic/examples/CTdataset/coords_lm.csv'
#serverPATH='./deepmedic/examples/CTdataset/Dataset/'
year='2012'
serverPATH='/media/ebermejo/My Book/archive/'+year

#2007 has 0 labeled sets
#2008 has 0 labeled sets
#2009 has 3 labeled sets
#2010 has 287 labeled sets:
#2011 has 3 labeled sets:
#2012 has 5 labeled sets:
continueSession=False
verbose=True

#Remove spanish accents
def normalize_character(s):
    s=unicode(s, errors='ignore')
    s2 = ''.join((c for c in unicodedata.normalize('NFD',unicode(s)) if unicodedata.category(c) != 'Mn'))
    return s2.decode()

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

#################### 3D landmark parsing method ####################
def parseLandmarks(filename):
    file=open(filename)
    count=0
    landmarkDict= collections.defaultdict(dict)
    landmarkList= []
    all=string.maketrans('','')
    nodigs=all.translate(all, string.digits)
    for line in file:
        count=count+1
        line=line.replace('"','')
        tokens=line.strip('\n').split(',')
        name=tokens[0].strip('\t ^') #.translate(all, nodigs) #In case we want to remove letters (name initials of patients)
        name=name.strip()
        if count==1:
            #Parse Landmark names
            for i in range(7,64,3):
                tk=tokens[i].split('.')
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
def readLastSession():
    series=[]
    paths=[]
    sfile=open('./last_series'+year+'.txt','r')
    pfile=open('./last_paths'+year+'.txt','r')
    scount=0
    pcount=0
    for line in sfile:
        scount=scount+1
        tokens=line.strip('\n').split(';')
        series.append([tokens[2]])
    for line in pfile:
        pcount=pcount+1
        tokens=line.strip('\n').split(';')
        paths.append(tokens[1]+'/')
    sfile.close()
    pfile.close()
    return series,paths

def searchPath(reader,input_dir):
    series=[]
    paths=[]
    tags=[]
    nbdl=1
    for cur_dir, subdirs, files in os.walk( input_dir ):
        if len(subdirs)==0:
            print('Read',cur_dir,'...')
            series_found = reader.GetGDCMSeriesIDs( cur_dir )
            if series_found:
                series.append( series_found )
                paths.append( cur_dir )
                tag = '_'.join( cur_dir.split( os.path.sep )[ nbdl-1 : ] ).replace(' ', '_')
                tags.append( tag )
    return series,paths,tags

if __name__ == '__main__':
    ################### Parse CSV file with 3D landmarks ####################
    land3D,landmarkList=parseLandmarks(csvFILE)
    print('Reading 3D landmark coordinates',csvFILE,'...')

    patients=Set([])
    ################### Initialize reading directories  ####################
    reader = sitk.ImageSeriesReader()
    if(continueSession):
        series,paths=readLastSession() #If the path was already processed
    else:
        series,paths,tags=searchPath(reader,serverPATH)# Begin searching the database structure and saving all the paths including a Dicom serie
        #The database is structured by year of acquisition (2007-2012), so we split the processing to avoid overheating the external drive
        sfile=open('./last_series'+year+'.txt','w')
        pfile=open('./last_paths'+year+'.txt','w')
        for ii, path in enumerate( paths ):
            pfile.write(str(ii)+';'+str(path)+'\n')
            for jj, serie in enumerate(series[ii]):
                sfile.write(str(ii)+';'+str(jj)+';'+str(serie)+'\n')
        pfile.close()
        sfile.close()
    ################### Process series: read image, obtain patient name and id, export dcm image to our own Database folder if labelled ###################
    for ii, path in enumerate( paths ):
        for jj, serie in enumerate(series[ii]):
            if(verbose):
                print (path,serie)
            try:
                dicom_names = reader.GetGDCMSeriesFileNames( path, serie)
            except:
                continue
            for zz in range(0,len(dicom_names)):
                exitFlag=False
                patientdir=patientname=patientid=''
                dcmf = dicom.read_file(dicom_names[zz])
                #Extract both Patient's ID and Name, as they identify our labelled studies
                patientid=dcmf.data_element('PatientID').value
                patientname=dcmf.data_element('PatientsName').value
                patientid=str(patientid).strip()
                patientname=str(patientname).strip().strip('\t ^')
                #Check if the patient is labelled
                if(land3D[patientname]):
                    patientdir=patientname
                elif(land3D[patientid]):
                    patientdir=patientid
                else:
                    print('Patient not recognized',patientid,patientname)
                    patients.add((dicom_names[zz],patientid,patientname,0)) #Save the path and the id
                    break
                #The patient has been recognized in our labelled dataset
                patients.add((dicom_names[zz],patientid,patientname,1)) #Save the path and the id for record
                #Extract series UID, description or study description to mimick Osirix export tool in the structure of our data folder
                serie=dcmf.data_element('SeriesInstanceUID').value
                try:
                    seriename=dcmf.data_element('SeriesDescription').value
                    seriename=seriename.replace(' ','').replace('.','')
                except:
                    try:
                        seriename=dcmf.data_element('StudyDescription').value
                        seriename=seriename.replace(' ','_').replace('.','')
                    except:
                        seriename='unknown'
                seriesnumber=dcmf.data_element('SeriesNumber').value
                sn=str(seriesnumber)
                seriename=normalize_character(seriename) #In case there are spanish accents in the descriptions
                seriedir=str(seriename)+'_'+sn
                imagedir='IM-'+sn.zfill(4)+'-'+'{:04d}'.format(zz+1)+'.dcm'
                if(verbose):
                    print(serie,patientdir,seriedir)
                #saving image in a predefined path style mimicking Osirix
                newdir='/media/ebermejo/My Passport/Dataset/'+year+'sets/'+patientdir+'/'+seriedir+'/'
                newfiledir=newdir+imagedir
                #Check if file was already exported [patient collisions when multiple studies involved]
                if not os.path.exists(newfiledir):
                    if not os.path.exists(newdir):
                        os.makedirs(newdir)
                    try:
                        dcmf.save_as(newfiledir)
                    except:
                        print('Check disk space, probably full',dicom_names[zz],newfiledir)
                        pause()
                    print('Saving',newfiledir,'...')
                else:
                    print('FILE exists already...skipping','...')
    ##Save the record of patient identifier and paths for future works
    namesfile=open('./last_patients'+year+'.txt','w')
    for ip in enumerate(patients):
        namesfile.write(str(ip)+'\n')
    namesfile.close()
