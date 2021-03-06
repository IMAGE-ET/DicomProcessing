import SimpleITK as sitk
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import skimage.io as io
import sys,os,fnmatch,shutil,unicodedata
import collections,random,operator
#################### Global constants ####################
verbose=True #used to control output progress information
verboseIMG=False #used to control visualization output
random.seed(14259)

#Image thresholds to purge undesired dcm series
min_landmarks_per_patient=1
min_size = 25 # voxels
max_voxel_spacing = 3.1 # mm
min_vert_extent = 10 #mm
#Auxiliary filenames for diferent paths and configuration files of Deepmedic
csvFILE='./coords_lm.csv'
datasetDIR='/pathToDataset/'
outputDIR = './processedData/'


#################### Managment tools ####################
def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def normalize_character(s):
    s=unicode(s, errors='ignore')
    s2 = ''.join((c for c in unicodedata.normalize('NFD',unicode(s)) if unicodedata.category(c) != 'Mn'))
    return s2.decode()

#################### Folder cleaning method ####################
def cleanTrainingFolder():
    if(os.path.exists(os.path.dirname(outputDIR))):
        while True:
            try:
                print('To avoid training data duplicates, train/validation/test folders must be empty')
                delete = raw_input('Should I remove the contents of these folders (y/n): ')
                #Python 3 -> change for input()
                #delete=str(delete)
                if(delete!='y' and delete!='n'):
                    raise ValueError
            except ValueError:
                print("Please enter y or n without quotes.")
                #better try again... Return to the start of the loop
                continue
            else:
                #Input successfully parsed!
                if(delete=='y'):
                    if os.path.exists(os.path.dirname(outputDIR)):
                        shutil.rmtree(outputDIR)
                break


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

#################### Show image intensity distribution ####################
def displayIntensity(image):
    image_size=image.GetSize()
    mask_ranges = [range(0,image_size[0]), range(0,image_size[1]), range(0, image_size[2])]
    intensity_values = sitk.GetArrayFromImage(image)
    roi_intensity_values = intensity_values[mask_ranges[2][0]:mask_ranges[2][-1],
                                          mask_ranges[1][0]:mask_ranges[1][-1],
                                          mask_ranges[0][0]:mask_ranges[0][-1]].flatten()
    plt.hist(roi_intensity_values, bins=100)
    plt.title("Intensity Values")
    plt.show()

######################   Preprocessing methods   ###############################
################################################################################
################################################################################
################################################################################

#################### Create Ground-Truth segmentation image ####################
def extractGT(image,roilabel,filename,land3D,landmarkList):
    ### Create a blank label image ###
    gt_label = sitk.Image(image.GetSize(), sitk.sitkUInt16)
    gt_label.SetOrigin(image.GetOrigin())
    gt_label.SetSpacing(image.GetSpacing())
    gt_label.SetDirection(image.GetDirection())
    nland=0
    count=0
    listi=[-1,0,1]
    listk=[0] #Draws a 3x3 square instead of a 3x3x3 cube using [-1,0,1]
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
                if(index3D[0]+i>=0 and index3D[0]+i<image.GetSize()[0]):
                    for j in listi:
                        if(index3D[1]+j>=0 and index3D[1]+j<image.GetSize()[1]):
                            for k in listk:
                                if(index3D[2]+k>=0 and index3D[2]+k<image.GetSize()[2]):
                                    #draws a "sphere" instead of a 3x3x3 box
                                    # if(k!=0 and i!=0 and j!=0):
                                    #     pass
                                    # else:
                                    block3D=[index3D[0]+i,index3D[1]+j,index3D[2]+k]
                                    if(roilabel[block3D]==1):
                                        gt_label[block3D]=it+1
            nland=nland+1
    if(verboseIMG):
        st='gt_overlay_'+str(nland)+'_'+dcmfilename
        img_T1_255 = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
        sitk.Show(sitk.LabelOverlay(img_T1_255,gt_label),st)
        pause()
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

def padImage(image,size):
    #Padding image centering it to facilitate posterior stiding using a 25x25x25 window
    bound_x=25-size[0]%25
    bound_y=25-size[1]%25
    bound_z=25-size[2]%25
    #Recalculate lower and upper bounds
    lb_x=int(bound_x/2.0)
    ub_x=int(bound_x-lb_x)
    lb_y=int(bound_y/2.0)
    ub_y=int(bound_y-lb_y)
    lb_z=int(bound_z/2.0)
    ub_z=int(bound_z-lb_z)
    #Run the Padding filter
    padder= sitk.ConstantPadImageFilter()
    padder.SetConstant(-2000)#As we are clamping and filtering the image to -2000, we doesn't need to get the minimum value of the image as background
    padder.SetPadUpperBound([ub_x,ub_y,ub_z])
    padder.SetPadLowerBound([lb_x,lb_y,lb_z])
    padded_image=padder.Execute(image)
    return padded_image

#################### Create image ROI from resampled image  ####################
def extractROI(image):
    ###  First we rescale the image for visualization     ###
    img_T1_255 = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
    size=img_T1_255.GetSize()
    ### And segment the background creating the label void ###
    thresholder=sitk.BinaryThresholdImageFilter()
    thresholder.SetInsideValue(0)
    thresholder.SetOutsideValue(1)
    thresholder.SetLowerThreshold(-3500) #Void (Black) intensity of the outer boundary (pan)
    #thresholder.SetUpperThreshold(-500) #Background intensity of the scans

    thresholder.SetUpperThreshold(-1000) #Air
    void_label=thresholder.Execute(image)

    ### Then we segment the brain white matter creating its own label ###
    thresholder.SetInsideValue(1)
    thresholder.SetOutsideValue(0)
    thresholder.SetLowerThreshold(20)
    thresholder.SetUpperThreshold(30)
    white_label=thresholder.Execute(image)

    ### We also segment the brain grey matter creating its own label ###
    thresholder.SetInsideValue(1)
    thresholder.SetOutsideValue(0)
    thresholder.SetLowerThreshold(37)
    thresholder.SetUpperThreshold(45)
    grey_label=thresholder.Execute(image)
    brain_label=white_label+grey_label

    ### Extract the image mask containing only our interest structures (bone+soft tissue) ###
    subject = sitk.BinaryMorphologicalClosing( void_label, 1 )
    brain_mask = sitk.BinaryMorphologicalClosing(brain_label, 1 )
    roi_label=subject#-brain_mask
    roi_mask = sitk.BinaryMorphologicalClosing( roi_label, 3 )
    roi_image=sitk.Mask(image,roi_mask)

    ### Show final segmentation overlay and slices ###
    if(verboseIMG):
        size=roi_image.GetSize()
        st='roi_mask'+'_'+dcmfilename
        #myshow3d(sitk.LabelOverlay(img_T1_255,roi_mask), zslices=range(50,size[2]-50,15), dpi=90, title=st)
        sitk.Show(sitk.LabelOverlay(img_T1_255,roi_mask),st)
        #sitk.Show(roi_mask)
    return roi_image,roi_mask

#################### Normalize Intensity Distribution ####################
def normalizeImg(image,roi):
    #applies pixel-wise a linear transformation to the intensity values of input image pixels.
    #The linear transformation is defined by the user in terms of the minimum and maximum values that the output image should have.
    #### TODO: Find proper intensity range?
    image_filter=sitk.ClampImageFilter()
    image_filter.SetLowerBound(-2000)
    image_filter.SetUpperBound(2000)
    filtered_image=image_filter.Execute(image)

    ### Using NormalizeImageFilter ###
    #shifts and scales an image so that the pixels in the image have a zero mean and unit variance.
    #This filter uses StatisticsImageFilter to compute the mean and variance of the input and
    #then applies ShiftScaleImageFilter to shift and scale the pixels.

    ### Using directly the ShiftScaleImageFilter ###
    #As we already have a known Housenfield distribution in our data,
    #apply Shift and Scale over our thresholded data
    shift_filter=sitk.ShiftScaleImageFilter()
    shift_filter.SetScale(1/4000.0)
    shift_filter.SetShift(2000)
    normalized_image=shift_filter.Execute(filtered_image)
    normalized_masked_image=sitk.Mask(normalized_image,roi)

    if(verboseIMG):
        sitk.Show(normalized_masked_image,'normalized_ct_image')
        displayIntensity(normalized_image)
        pause()
    return normalized_image,normalized_masked_image
#################### Bone structure segmentation ####################
#### NOT USED FOR NOW ###
# def extractBone(image):
#     ###  First we rescale the image for visualization     ###
#     img_T1_255 = sitk.Cast(sitk.RescaleIntensity(image), sitk.sitkUInt8)
#     size=img_T1_255.GetSize()
#     thresholder=sitk.BinaryThresholdImageFilter()
#     thresholder.SetInsideValue(1)
#     thresholder.SetOutsideValue(0)
#     thresholder.SetLowerThreshold(400)
#     thresholder.SetUpperThreshold(4000)
#     bone_label=thresholder.Execute(image)
#     bone_label = sitk.BinaryMorphologicalClosing( bone_label, 1 )
#     ### image containing only bone structure ###
#     bone_image=sitk.Mask(image,bone_label)
#     ### Show bone segmentation overlay and slices ###
#     if(verboseIMG):
#         size=bone_label.GetSize()
#         myshow3d(sitk.LabelOverlay(img_T1_255,bone_label), yslices=range(50,size[1]-50,15), zslices=range(50,size[2]-50,15), dpi=90,title='bone_slices')
#         sitk.Show(bone_image,'bone_image')
#     return bone_image,bone_label

#################### Image processing steps ####################
def preprocessImage(input_image,itcounter,land3D,landmarkList):
    #################### Get useful image data ####################
    spacing = input_image.GetSpacing()
    size = input_image.GetSize()
    origin = input_image.GetOrigin()
    direction = input_image.GetDirection()
    extent = [ size[j]*spacing[j] for j in range(3) ]

    if(verbose):
        print('Input_size:',size)
        print('Input_spacing',spacing)
        print('Input_extent',extent)
    #################### Recalculate resampled image size and spacing ####################
    if(verbose):
        print('Resampling to Isotropic image','...')
    new_size=[int(size[0]*(spacing[0])),int(size[1]*(spacing[1])),int(size[2]*(spacing[2]))]
    new_spacing = [1.0,1.0,1.0]
    isotropic_image=resampleImage(input_image,new_size,origin,direction,new_spacing)
    padded_image=padImage(isotropic_image,new_size)
    #################### Create image ROI from resampled image  ####################
    if(verbose):
        print('Extracting ROI','...')
    roi_image,roi_mask=extractROI(padded_image)
    #bone_image,bone_mask=extractBone(isotropic_image)

    #################### Create Ground-Truth segmentation ####################
    if(verbose):
        print('Generating Ground-Truth labels','...')
    ##Test image with original size -> same results
    gt_image=extractGT(padded_image,roi_mask,dcmfilename,land3D,landmarkList)

    #################### Normalize Intensity Distribution ####################
    if(verbose):
        print('Normalizing Intensity Distribution','...')
    ### Rescale intensity distribution ###
    normalized_image,normalized_masked_image=normalizeImg(padded_image,roi_mask)
#################### Write NiFTi images #####################

    name_file_ct=dcmfilename+'/'+dcmfilename+'_ct.nii.gz'
    name_file_gt=dcmfilename+'/'+dcmfilename+'_gt.nii.gz'
    name_file_roi=dcmfilename+'/'+dcmfilename+'_roi.nii.gz'
    if(verbose):
        print('Writing patient "'+dcmfilename+'" NiFTi images','...')

    if not os.path.exists(os.path.dirname(outputDIR+name_file_ct)):
        try:
            os.makedirs(os.path.dirname(outputDIR+name_file_ct))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    sitk.WriteImage(normalized_masked_image,outputDIR+name_file_ct)
    sitk.WriteImage(roi_mask,outputDIR+name_file_roi)
    sitk.WriteImage(gt_image,outputDIR+name_file_gt)
    return 1,normalized_image.GetSize()
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':
    ### Ask user if train/validation/test folders are to be removed ###
    cleanTrainingFolder()

    log_file = open('processed_files.log', "w")
    #################### Parse CSV file with 3D landmarks ####################
    land3D,landmarkList,patientData=parseLandmarks(csvFILE)
    print('Reading 3D landmark coordinates',csvFILE,'...')

    #################### Parse Dicom Dataset directory    ####################
    ### unique dcm folders ###
    subjects=os.walk(datasetDIR).next()[1]
    ### Randomize the database for file classification in train,test or validation sets ###
    #nSubjects=len(land3D) #guided by the .csv files is safer (as some images have no landmarks)
    nSubjects=len(subjects) #guided by the number of images

    reader = sitk.ImageSeriesReader()
    counter=-1
    preprocessed_images=[]
    preprocessed_sizes=[]
    for itsubject,subject in enumerate(subjects):
        scans,directories=searchPath(datasetDIR+subject)
        print('Reading',len(directories),'Dicom images',directories,'...')
        idealvolume=''
        optionalvolumes=[]
        optionalsizes=[]
        ### Iterate subject directory series in order to find the suitable volume ###
        for itfile,dirname in enumerate(directories):
            string=dirname.replace(datasetDIR,'')
            string=string.strip('.').split('/')
            dcmfilename=string[0]
            dcmserie=string[1]
            if(len(string)>2):
                dcmserie=dcmserie+'_'+string[2]
            dicom_names = reader.GetGDCMSeriesFileNames( dirname )
            reader.SetFileNames(dicom_names)
            try:
                test_image = reader.Execute()
            except:
                print('Failed to read image, please check',dcmfilename,dcmserie,'...')
                continue
            #################### Get useful image data ####################
            spacing = test_image.GetSpacing()
            size = test_image.GetSize()
            origin = test_image.GetOrigin()
            direction = test_image.GetDirection()
            extent = [ size[j]*spacing[j] for j in range(3) ]

            if min( size ) < min_size:
                print('Skipping due to size',itfile,dcmfilename,dcmserie,size,'...')
                continue
            if max( spacing ) > max_voxel_spacing:
                print('Skipping due to spacing',itfile,dcmfilename,dcmserie,spacing,'...')
                continue
            if extent[2] < min_vert_extent:
                print('Skipping due to extent',itfile,dcmfilename,dcmserie,extent,'...')
                continue
            optionalvolumes.append(dirname)
            optionalsizes.append(size[2])

        if(not optionalvolumes):
            continue
        index, value = max(enumerate(optionalsizes), key=operator.itemgetter(1))

        if(optionalvolumes[index]):
            idealvolume=optionalvolumes[index]

        ### Once we found a suitable volume, process it ###
        if(idealvolume):
            string=idealvolume.replace(datasetDIR,'')
            string=string.strip('.').split('/')
            dcmfilename=string[0]
            dcmserie=string[1]
            if(len(string)>2):
                dcmserie=dcmserie+'_'+string[2]
            has3Dlandmark=False

            if(land3D[dcmfilename]):
                for it,land in enumerate(landmarkList):
                    point3D=land3D[dcmfilename][landmarkList[it]]
                    lx,ly,lz=point3D
                    if(lx!='' and ly!='' and lz!=''):
                        has3Dlandmark=True
            if(has3Dlandmark):
                counter=counter+1
                dicom_names = reader.GetGDCMSeriesFileNames( idealvolume )
                reader.SetFileNames(dicom_names)
                input_image = reader.Execute()
                print('Preprocessing',counter,dcmfilename,dcmserie,'...')
                result,output_size=preprocessImage(input_image,counter,land3D,landmarkList)

                if( result==1):
                    print('Preprocessing successful',counter,'...')
                    preprocessed_images.append(idealvolume)
                    preprocessed_sizes.append(output_size)
                else:
                    print('Something happened','...')
            else:
                print('Subject',dcmfilename,'was skipped due not having any 3D landmark','...')
        else:
            print(dcmfilename,' could not be processed, check image properties','...')

    print('Database preprocessing successful. Number of images',counter,'items','...')
    s="{s}-[{d1},{d2},{d3}]\n"
    for index,item in enumerate(preprocessed_images):
        log_file.write(s.format(s=item,d1=preprocessed_sizes[index][0],d2=preprocessed_sizes[index][1],d3=preprocessed_sizes[index][2]))
    log_file.close()


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
