
import SimpleITK as sitk
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import skimage.io as io
import sys,os,fnmatch,shutil,string,unicodedata
import collections,random
import dicom
from nibabel.testing import data_path
from ipywidgets import interact,fixed
from distutils.dir_util import copy_tree
#################### Global constant ####################
#Relative paths to our dataset and output files
year='2010'
#datasetPATH='./Dataset/'
datasetPATH='/media/ebermejo/My Passport/Dataset/'+year+'sets/'
csvFILE='./deepmedic/examples/CTdataset/coords_lm.csv'
outputPATH='./cleanDataset/'
log_file = open('processed_files'+year+'.log', "w")

log_file_disc = open('processed_files'+year+'_discarded.log', "w")

#Image thresholds to purge undesired dcm series
min_landmarks_per_patient=1
min_size = 5 # voxels
max_voxel_spacing = 5.1 # mm
#min_vert_extent = 50.0 # mm

#Other useful variables
yes=set(['y','','yes','ye'])
no=set(['n','no'])
verbose=True #used to control output progress information
verboseIMG=True #used to control visualization output

#################### Managment tools ####################
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def normalize_character(s):
    s=unicode(s, errors='ignore')
    s2 = ''.join((c for c in unicodedata.normalize('NFD',unicode(s)) if unicodedata.category(c) != 'Mn'))
    return s2.decode()

#################### 3D landmark parsing method for our labelled data in a csv file ####################
#csv file includes a patient identifier, followed by age, age group and some data relative to the type of
#the scan. Then 3D coordinates are included for each of the 19 anatomical landmarks considered in our study
def parseLandmarks(filename):
    file=open(filename)
    count=0
    landmarkDict= collections.defaultdict(dict)
    landmarkList= []
    patientData= collections.defaultdict(dict)
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
            data=normalize_character(tokens[6].strip())
            patientData[name]=[tokens[1].zfill(3),str(data)]
            countl=0
            for i in range(7,64,3):
                lx=tokens[i].strip()
                ly=tokens[i+1].strip()
                lz=tokens[i+2].strip()
                if(lx!='' and ly!='' and lz!=''):
                    countl+=1
            #use the previous loop in case we dont want to save blank landmarks in the dictionary
            if(countl>=min_landmarks_per_patient):
                for i in range(7,64,3):
                    lx=tokens[i]
                    ly=tokens[i+1]
                    lz=tokens[i+2]
                    try:
                        landmarkDict[name][landmarkList[listindex]]=[lx,ly,lz]
                    except:
                        print(name,'ERROR adding name')
                    listindex=listindex+1
    print('Number of 3D Landmarks',len(landmarkList))
    print('Number of Subjects',len(landmarkDict),len(patientData))
    return landmarkDict,landmarkList,patientData
    file.close()

#################### Show methods definition for plotting slices ####################
def myshow(img, title=None, margin=0.05, dpi=80 ):
    nda = sitk.GetArrayFromImage(img)
    spacing = img.GetSpacing()
    if nda.ndim == 3:
        c = nda.shape[-1]
        ### if the number of components is 3 or 4 consider it an RGB image ###
        if not c in (3,4):
            nda = nda[nda.shape[0]//2,:,:]
    elif nda.ndim == 4:
        c = nda.shape[-1]
        if not c in (3,4):
            raise Runtime("Unable to show 3D-vector Image")
        ### take a z-slice ###
        nda = nda[nda.shape[0]//2,:,:,:]
    ysize = nda.shape[0]
    xsize = nda.shape[1]
    ### Make a figure big enough to accomodate an axis of xpixels by ypixels ###
    ### as well as the ticklabels, etc. ###
    figsize = (1 + margin) * ysize / dpi, (1 + margin) * xsize / dpi
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ### Make the axis the right size ###
    ax = fig.add_axes([margin, margin, 1 - 2*margin, 1 - 2*margin])
    extent = (0, xsize*spacing[1], ysize*spacing[0], 0)
    t = ax.imshow(nda,extent=extent,interpolation=None)
    if nda.ndim == 2:
        t.set_cmap("gray")
    if title:
        plt.title(title)
    plt.show()

def myshow3d(img, xslices=[], yslices=[], zslices=[], title=None, margin=0.05, dpi=80):
    size = img.GetSize()
    img_xslices = [img[s,:,:] for s in xslices]
    img_yslices = [img[:,s,:] for s in yslices]
    img_zslices = [img[:,:,s] for s in zslices]
    maxlen = max(len(img_xslices), len(img_yslices), len(img_zslices))
    img_null = sitk.Image([0,0], img.GetPixelIDValue(), img.GetNumberOfComponentsPerPixel())
    img_slices = []
    d = 0
    if len(img_xslices):
        img_slices += img_xslices + [img_null]*(maxlen-len(img_xslices))
        d += 1
    if len(img_yslices):
        img_slices += img_yslices + [img_null]*(maxlen-len(img_yslices))
        d += 1
    if len(img_zslices):
        img_slices += img_zslices + [img_null]*(maxlen-len(img_zslices))
        d +=1
    if maxlen != 0:
        if img.GetNumberOfComponentsPerPixel() == 1:
            img = sitk.Tile(img_slices, [maxlen,d])
        else:
            img_comps = []
            for i in range(0,img.GetNumberOfComponentsPerPixel()):
                img_slices_c = [sitk.VectorIndexSelectionCast(s, i) for s in img_slices]
                img_comps.append(sitk.Tile(img_slices_c, [maxlen,d]))
            img = sitk.Compose(img_comps)
    myshow(img, title, margin, dpi)

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

######################   Preprocessing methods   ###############################
################################################################################
################################################################################
################################################################################

#################### Create Ground-Truth segmentation image ####################
def extractGT(image,filename,land3D,landmarkList):
    ### Create a blank label image ###
    gt_label = sitk.Image(image.GetSize(), sitk.sitkUInt16)
    gt_label.SetOrigin(image.GetOrigin())
    gt_label.SetSpacing(image.GetSpacing())
    gt_label.SetDirection(image.GetDirection())
    nland=0
    count=0
    listi=[-1,0,1]
    ### Iterate the list of landmarks ###
    for it,land in enumerate(landmarkList):
        point3D=land3D[filename][landmarkList[it]]
        lx,ly,lz=point3D
        if(lx!='' and ly!='' and lz!=''):
            ### Assign the landmark label to the GT image       ###
            ### Note that label index must start from 1 and     ###
            ### increasing by one as Background is labeled as 0 ###
            ### Besides we have to find the voxel of the        ###
            ### 3D phisical point from the csv file             ###
            point3D=[float(lx),float(ly),float(lz)]
            voxel3D=image.TransformPhysicalPointToIndex(point3D)
            index3D=[int(voxel3D[0]),int(voxel3D[1]),int(voxel3D[2])]
            count=count+1
            for i in listi:
                for j in listi:
                    for k in listi:
                        block3D=[index3D[0]+i,index3D[1]+j,index3D[2]+k]
                        try:
                            gt_label[block3D]=it+1
                        except:
                            print('Skipping serie due to different area in scan',filename,'...')
                            return 0
            nland=nland+1
    #Visualize the labeled groundtruth to confirm it is correct
    if(verboseIMG):
        st='gt_overlay_'+str(nland)+'_'+dcmfilename
        img_T1_255 = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
        sitk.Show(sitk.LabelOverlay(img_T1_255,gt_label),st)
    return gt_label

################# Create resampler and obtain isotropic_image  #################
def resampleImage(image,new_size,origin,direction,new_spacing):
    resampler = sitk.ResampleImageFilter()
    resampler.SetSize(new_size)
    resampler.SetOutputOrigin(origin)
    resampler.SetOutputDirection(direction)
    resampler.SetOutputSpacing(new_spacing)
    resampled_image = resampler.Execute(image)
    return resampled_image


#################### Image processing steps ####################
def preprocessImage(input_image,land3D,landmarkList,dirname,dcmfilename,dcmserie):
    global counter
    saved=-1
    #################### Get useful image data ####################
    spacing = input_image.GetSpacing()
    size = input_image.GetSize()
    origin = input_image.GetOrigin()
    direction = input_image.GetDirection()
    extent = [ size[j]*spacing[j] for j in range(3) ]
    #print('Pixel0',image.GetPixel(0, 0, 0))
    if(verbose):
        print('Input_size:',size)
        print('Input_spacing',spacing)
        print('Input_extent',extent)
    #################### Recalculate resampled image size and spacing ####################
    # if(verbose):
    #     print('Resampling to Isotropic image','...')
    # new_size=[int(size[0]*(spacing[0])),int(size[1]*(spacing[1])),int(size[2]*(spacing[2]))]
    # new_spacing = [1.0,1.0,1.0]
    # isotropic_image=resampleImage(input_image,new_size,origin,direction,new_spacing)

    #################### Create Ground-Truth segmentation ####################
    if(verbose):
        print('Generating Ground-Truth labels','...')
    ##Test image with original size -> same results
    gt_image=extractGT(input_image,dcmfilename,land3D,landmarkList)
    if(gt_image):
        pause()
        print('Please, verify the landmarks are placed porperly.','...')
        print('Validating serie from',dirname,'...')
        done=False
        while(done==False):
            check=raw_input(' Do you want to validate this serie (y/n):')
            if (check in yes):
                done=True
                counter+=1
                sourceDCMPath=dirname
                outputDCMPath=outputPATH+dcmfilename+'/'+dcmserie+'/'
                print('Saving to:',outputDCMPath,'...')
                # if not os.path.exists(os.path.dirname(outputDCMPath)):
                #     try:
                #         os.makedirs(os.path.dirname(outputDCMPath))
                #     except OSError as exc: # Guard against race condition
                #         if exc.errno != errno.EEXIST:
                #             raise
                copy_tree(sourceDCMPath, outputDCMPath,update=0,verbose=1)
                saved=2
            elif(check in no):
                done=True
                saved=1
            else:
                done=False
                sys.stdout.write('Please use y/yes/intro or no/n')
    else:
        saved=1
    return saved
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':
    global counter
    counter=0
    #################### Parse CSV file with 3D landmarks ####################
    land3D,landmarkList,patientData=parseLandmarks(csvFILE)
    print('Reading 3D landmark coordinates',csvFILE,'...')

    #################### Parse Dicom Dataset directory    ####################
    ### unique dcm folders ###
    subjects=os.walk(datasetPATH).next()[1]

    reader = sitk.ImageSeriesReader()
    preprocessed_images=[]
    discarded_images=[]
    for itsubject,subject in enumerate(subjects):
        scans,directories=searchPath(datasetPATH+subject)
        print('Reading',len(directories),'Dicom images',directories,'...')
        optionalvolumes=[]
    #    forgetSubject=False
        #Verifiy and count the landmarks of this subject
        has3Dlandmark=False
        countl=0
        if(land3D[subject]):
            for it,land in enumerate(landmarkList):
                    point3D=land3D[subject][landmarkList[it]]
                    lx,ly,lz=point3D
                    if(lx!='' and ly!='' and lz!=''):
                        has3Dlandmark=True
                        countl+=1
        else:
            print('Subject',subject,'was skipped due to not having any 3D landmark','...')
        if(has3Dlandmark):
            ### Iterate subject directory series in order to find a subset of suitable volumes ###
            for itfile,dirname in enumerate(directories):
                #Extract patientID and serieDescription from directory
                string=dirname.replace(datasetPATH,'')
                string=string.strip('.').split('/')
                dcmfilename=string[0] #<-aka subject
                dcmserie=string[1]
                if(len(string)>2):
                    dcmserie=dcmserie+'_'+string[2]
                #First of all, check the Age of the subject in case we have duplicates.
                try:
                    _,_,filenames=os.walk(dirname).next()
                    dcmf = dicom.read_file(dirname+'/'+filenames[0],force=True)
                    patientage=dcmf.data_element('PatientsAge').value.strip('Y')
                    if(patientage!=patientData[dcmfilename][0]):
                        print("Anomaly found, patient age does not correspond with label, please verify")
                        print(subject,patientage,patientData[dcmfilename][0])
                        discarded_images.append([dcmfilename,dcmserie])
                        continue
                except:
                    print('Age could not be confirmed','...')
                #Skip series known to be crops of original scan
                if(('Axial' in dcmserie) or ('Coronal' in dcmserie) or ('Sagittal' in dcmserie)):
                    if(verbose):
                        print('Skipping due to Axial/Coronal/Saggital crop',itfile,dcmfilename,dcmserie,size,'...')
                        discarded_images.append([dcmfilename,dcmserie])
                    continue

                #Read the image to extract data
                dicom_names = reader.GetGDCMSeriesFileNames( dirname )
                reader.SetFileNames(dicom_names)
                try:
                    test_image = reader.Execute()
                except:
                    print('Failed to read image, please check before next attempt',dcmfilename,dcmserie,'...')
                    pause()
                    continue

                #################### Get useful image data ####################
                spacing = test_image.GetSpacing()
                size = test_image.GetSize()
                extent = [ size[j]*spacing[j] for j in range(3) ]

                #Skip images with small number of slices
                if min( size ) < min_size:
                    print('Skipping due to size',itfile,dcmfilename,dcmserie,size,'...')
                    discarded_images.append([dcmfilename,dcmserie])
                    continue
                #Skip images with spacing too wide
                if max( spacing ) > max_voxel_spacing:
                    print('Skipping due to spacing',itfile,dcmfilename,dcmserie,spacing,'...')
                    discarded_images.append([dcmfilename,dcmserie])
                    continue
                # if extent[2] < min_vert_extent:
                #     if(verbose):
                #         print('Skipping due to extent',itfile,dcmfilename,dcmserie,extent,'...')
                #     continue
                print('Added series for deep processing',dirname,'...')
                optionalvolumes.append(dirname)
            if(optionalvolumes):
                for itfile,dirname in enumerate(optionalvolumes):
                    string=dirname.replace(datasetPATH,'')
                    string=string.strip('.').split('/')
                    dcmfilename=string[0] #<-aka subject
                    dcmserie=string[1]
                    if(len(string)>2):
                        dcmserie=dcmserie+'_'+string[2]
                    dicom_names = reader.GetGDCMSeriesFileNames( dirname )
                    reader.SetFileNames(dicom_names)
                    try:
                        test_image = reader.Execute()
                    except:
                        continue
                    #################### Get useful image data ####################
                    spacing = test_image.GetSpacing()
                    size = test_image.GetSize()

                    print('Subject has',countl,'landmarks','...')
                    dicom_names = reader.GetGDCMSeriesFileNames( dirname )
                    reader.SetFileNames(dicom_names)
                    input_image = reader.Execute()
                    print('Preprocessing',dcmfilename,dcmserie,'...')
                    result=preprocessImage(input_image,land3D,landmarkList,dirname,dcmfilename,dcmserie)

                    if( result==2):
                        print('Preprocessing successful','...')
                        preprocessed_images.append([dcmfilename,dcmserie])
                    elif(result==1):
                        print('Image not saved','...')
                        discarded_images.append([dcmfilename,dcmserie])
                    else:
                        print('Something happened during preprocessing','...')
                        pause()

            else:
                print(dcmfilename,' could not be processed, check image properties','...')
            # programPause = raw_input("Press the <ENTER> key to continue...")

    print('Database preprocessing successful. Number of images',counter,'items','...')
    for item in preprocessed_images:
        log_file.write("%s\n" % item)
    log_file.close()
    for item in discarded_images:
        log_file_disc.write("%s\n"%item)
    log_file_disc.close()
#################### edge detection ####################
# edges = sitk.CannyEdgeDetection(sitk.Cast(isotropic_image, sitk.sitkFloat32), lowerThreshold=0.0,
#                                 upperThreshold=200.0, variance = (5.0,5.0,5.0))
# edge_indexes = np.where(sitk.GetArrayFromImage(edges) == 1.0)
#
# # Note the reversed order of access between SimpleITK and numpy (z,y,x)
# physical_points = [edges.TransformIndexToPhysicalPoint([int(x), int(y), int(z)]) \
#                    for z,y,x in zip(edge_indexes[0], edge_indexes[1], edge_indexes[2])]
#
# edge_label = sitk.Image(isotropic_image.GetSize(), sitk.sitkUInt16)
# edge_label.CopyInformation(isotropic_image)
# e_label = 255
# for point in physical_points:
#     edge_label[edge_label.TransformPhysicalPointToIndex(point)] = e_label
#
# sitk.Show(sitk.LabelOverlay(sitk.Cast(sitk.IntensityWindowing(isotropic_image, windowMinimum=-32767, windowMaximum=-29611),
#                                       sitk.sitkUInt8), edge_label, opacity=0.5))
