import SimpleITK as sitk
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import skimage.io as io
import sys,os,fnmatch,shutil
import collections,random
from nibabel.testing import data_path
from ipywidgets import interact,fixed
import re
import scipy.spatial as spatial


#################### Global variables ####################
verbose=False #used to control output progress information
verboseIMG=True #used to control visualization output

datasetDIR = './'
csvFILE='coords_lm.csv'
resultFILE='./SingleLand3D_results.csv'
distanceFILE='./SingleLand3D_distances.csv'
outputDIR = './'
# testDIR='./setEctos/test/'
# predictionsDIR='../resultsEctos/predictions/testSessionEctoLandmarks/predictions/'
testDIR='./set19Lands/test/'
predictionsDIR='../resultsSingleNoVa/predictions/testSessionSingle3DCnn/predictions/'
#predictionsDIR='../SingleLandoutput3D/predictions/testSessionSingle3DCnn/predictions/'
n_classes=19

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]
def centeroidnp(arr):
    print(arr.shape)
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return int(sum_x/length), int(sum_y/length), int(sum_z/length)

####################  Path reading methods for NifTi image predictionss associated with the tested image:  ####################
def searchPath(input_dir,fileid='*.nii.gz'):
    matches = []
    dirs = []
    for root, dirnames, filenames in os.walk(input_dir):
        for filename in fnmatch.filter(filenames, fileid):
            matches.append(os.path.join(root, filename));
            if root not in dirs:
                dirs.append(root)
    return (matches, dirs)

#################### Predicted landmarks extraction method ####################
def neighboring( array ):
    nn,mm = len(array), len(array[0])
    offset = (0,-1,1) # 0 first so the current cell is the first in the gen
    indices = ( (i,j) for i in range(nn) for j in range(mm) )
    for i,j in indices:
        all_neigh =  ( (i+x,j+y) for x in offset for y in offset )
        valid = ( (i,j) for i,j in all_neigh if (0<=i<nn) and (0<=j<mm) ) # -1 is a valid index in normal lists, but not here so throw it out
        yield valid.next(), valid ## first is the current cell, next are the neightbors


def extractPredictions(subject,image):
    image_array = sitk.GetArrayFromImage(image)
    size=np.shape(image_array)
    resultclass=[]
    resultclass=[(image_array[i][j][k],i,j,k) for i in range(size[0]) for j in range(size[1]) for k in range(size[1]) if image_array[i][j][k]!=0]
    sizerc=len(resultclass)
    listgt='gt_'+subject+','
    listpred='pred_'+subject+','
    listdist='dist_'+subject+','
    listdist2=''+subject+','
    for nclass in range(1,n_classes):
        dist=-1
        predcoord3D=np.array([-1,-1,-1])
        point3D=np.array(land3D[subject][landmarkList[nclass]])
        lx,ly,lz=point3D
        if(lx!='' and ly!='' and lz!=''):
            gtpoint3D=np.array([float(lx),float(ly),float(lz)])
        else:
            gtpoint3D=np.array([-1,-1,-1])
        resultclassi=[]
        resultclassi=[(resultclass[i][1],resultclass[i][2],resultclass[i][3]) for i in range(sizerc) if resultclass[i][0]==nclass]
        print(resultclassi)

        if resultclassi:
            a = np.array(resultclassi)
            print(a)
            point_tree = spatial.cKDTree(a)
            maxcluster=-1
            biggest=[]
            for center, group  in zip(a, point_tree.query_ball_point(a, 10)):
                cluster = point_tree.data[group]
                if(len(cluster)>maxcluster):
                    maxcluster=len(cluster)
                    biggest=center
            #if biggest:
            group= point_tree.query_ball_point(biggest, 10)
            cluster = point_tree.data[group]
            print(biggest,group)
            centeroid=centeroidnp(cluster)
            realcenter=[centeroid[2],centeroid[1],centeroid[0]]
            point3D=image.TransformIndexToPhysicalPoint(realcenter)
            lx,ly,lz=point3D
            predcoord3D=np.array([float(lx),float(ly),float(lz)])

            print(centeroid,predcoord3D,gtpoint3D)
            dist = np.linalg.norm(predcoord3D-gtpoint3D)
            print(predcoord3D,gtpoint3D,dist)
        if(predcoord3D[0]==-1 and predcoord3D[1]==-1 and predcoord3D[2]==-1):
            listpred+=',,,'
        else:
            listpred+=str(round(predcoord3D[0],1))+','+str(round(predcoord3D[1],1))+','+str(round(predcoord3D[2],1))+','
        if(gtpoint3D[0]==-1 and gtpoint3D[1]==-1 and gtpoint3D[2]==-1):
            listgt+=',,,'
        else:
            listgt+=str(round(gtpoint3D[0],1))+','+str(round(gtpoint3D[1],1))+','+str(round(gtpoint3D[2],1))+','
        if(dist==-1):
            listdist+=',,,'
            listdist2+=','
        else:
            listdist+=','+str(round(dist,2))+',,'
            listdist2+=str(round(dist,2))+','
    listgt+='\n'
    listpred+='\n'
    listdist+='\n'
    listdist2+='\n'
    outputFile.write(listgt)
    outputFile.write(listpred)
    outputFile.write(listdist)
    outputdistFile.write(listdist2)

#################### Predicted image selection method ####################
def visualizePredictions(subject,image,predictions):
    img_T1_255 = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
    for i,pred in enumerate(predictions):
        #idland=pred.replace(predictionsDIR,'').replace('pred_'+subject+'_ProbMapClass','').replace('.nii.gz','')
        if('Segm' in pred):
            reader=sitk.ImageFileReader()
            reader.SetFileName(pred)
            segmented_image = reader.Execute()
            extractPredictions(subject,segmented_image)
    #####  In case we prefer to stich each probability map altogether to show the output via probability map#####
    #     if('Segm' not in pred):
    #         reader.SetFileName(pred)
    #         imagepredicted= reader.Execute()
    #         if('Class0' in pred):
    #             st='Test_background_'+subject
    #             #sitk.Show(imagepredicted,st)
    #         elif('ProbMapClass' in pred):
    #             thresholder=sitk.BinaryThresholdImageFilter()
    #             thresholder.SetInsideValue(0)
    #             thresholder.SetOutsideValue(int(idland))
    #             thresholder.SetLowerThreshold(0) #Void (Black) intensity of the outer boundary (pan)
    #             thresholder.SetUpperThreshold(0.09) #Baclground intensity of the scans
    #             partial_label=thresholder.Execute(imagepredicted)
    #             if(idland==1 or i==0):
    #                 pred_label=partial_label
    #             else:
    #                 pred_label=pred_label+partial_label
    if(verboseIMG):
        st='Test_predictions_'+subject
        sitk.Show(sitk.LabelOverlay(img_T1_255,segmented_image),st)
        programPause = raw_input("Press the <ENTER> key to continue...")

#################### 3D landmark parsing method ####################
def parseLandmarks(filename):
    fileland=open(filename)
    count=0
    landmarkDict= collections.defaultdict(dict)
    landmarkList= []
    for line in fileland:
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
    fileland.close()

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
    listi=[-1,0,1]
    ### Iterate the list of landmarks ###
    for it,land in enumerate(landmarkList):
        if(1):
        #if(it==7): #We only pick one class this time using landmark Ec
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
                for i in listi:
                    for j in listi:
                        for k in listi:
                            block3D=[index3D[0]+i,index3D[1]+j,index3D[2]+k]
                            gt_label[block3D]=it+1
                nland=nland+1
    if(verboseIMG):
        st='gt_overlay_'+str(nland)+'_'+filename
        img_T1_255 = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
        sitk.Show(sitk.LabelOverlay(img_T1_255,gt_label),st)
    return gt_label


if __name__ == '__main__':
    #################### Parse CSV file with 3D landmarks ####################
    land3D,landmarkList=parseLandmarks(csvFILE)
    print('Reading 3D landmark coordinates',csvFILE,'...')

    #################### Parse Test Dataset directory    ####################
    ### subject id location ###
    subjects=os.walk(testDIR).next()[1]

    outputdistFile = open(distanceFILE, "w")
    outputFile = open(resultFILE, "w")
    header='Subject,'
    header2='Subject,'
    for i,land in enumerate(landmarkList):
        header+=land+'.x,'
        header+=land+'.y,'
        header+=land+'.z,'
        header2+=land+','
    header+='\n'
    header2+='\n'
    outputFile.write(header)
    outputdistFile.write(header2)

    reader = sitk.ImageSeriesReader() #Change to read nii images
    counter=-1
    preprocessed_images=[]
    for itsubject,subject in enumerate(subjects):
        scans,_=searchPath(testDIR+subject)
        for i,scan in enumerate(scans):
            if('ct' in scan):
                input_image=scan
            elif('gt in scan'):
                gt_image=scan
        print('Reading test image',input_image,'...')
        reader=sitk.ImageFileReader()
        reader.SetFileName(input_image)
        image = reader.Execute()
        predictions,_=searchPath(predictionsDIR,'*'+subject+'*')
        predictions.sort(key=natural_keys)
        ids=[]
        for x in predictions:
            ids.append(x.strip(predictionsDIR))
        print('Predictions,',ids)
        visualizePredictions(subject,image,predictions)
        extractGT(image,subject,land3D,landmarkList)
    outputFile.close()
