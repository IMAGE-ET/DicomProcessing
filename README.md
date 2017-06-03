# DicomProcessing
Dicom imaging and Deep learning processing set of tools for a specific database of CT scans used to extract anatomic landmarks
--

_______

Step 1 - Analize the server database path structure, compare the patients Id / name with our previously labelled information and extract the different series associated with each patient to a designated location

Step 2 - Semiautomatic refinement of the database for removing undesired series. Most of the series are automatically discarded, but visual confirmation for a proper landmark location is required.

Step 3 - Automatic preprocessing of the database:
        -Image resampling   - normalize image spacing to (1.0,1.0,1.0) mm per volex
        -ROI extraction     - creates a mask where blank space- air is not considered 
        -GT extraction      - creates a labelled mask where a 3x3x3voxels cube encloses a particular craniofacial landmark                                 (19 different classes)
        -Iage normalization - normalize intensity distribution for every image
        -NiFTi extraction   - saves the dicom images in NiFTi format using nibabel
        -DeepMedic configuration files creation
        
Step 4 - Visualization and error distances calculation for the predicted landmarks in the test dataset
