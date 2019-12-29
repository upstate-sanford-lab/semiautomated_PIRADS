import os
import re
import numpy as np
import pandas as pd
import pydicom
import shutil
import statistics
from collections import Counter
from datetime import date
import nibabel as nib

#author @t_sanf

class DatasetInterrogate:
    '''misc functions to help perform data analysis'''

    def __init__(self):
        self.basePATH='/home/tom/Desktop'
        self.database='revision_analysis2'
        self.devfolder='model_dev_indvPIRADS'

    def check_pt_spillage(self):
        '''function to ensure the training dataset did not spill over into the validation or testing dataset '''

        path=os.path.join(self.basePATH, self.database, self.devfolder)
        unique_mrns={}
        for db in ('train','val','test'):
            unique_mrns[db]=set(file.split('_')[0] for file in get_slices(path=path,dir=db))  #get_slices --> helper function at bottom
            print('{} unique mrns in database {}'.format(len(unique_mrns[db]),db))

        print("mrns overlap between training and validation sets: {}".format(unique_mrns['train'].intersection(unique_mrns['val'])))
        print("mrns overlap between training and test sets: {}".format(unique_mrns['train'].intersection(unique_mrns['test'])))
        print("mrns overlap between val and test sets: {}".format(unique_mrns['val'].intersection(unique_mrns['test'])))

    def ID_label_errors(self,ext='voi'):
        '''recursively looks through filetree and finds files that have not been labeled according to our conventions'''
        problems=[]
        for root, dirnames, filenames in os.walk(self.basePATH):
            for filename in filenames:
                if filename.endswith(ext):
                    filename_list=filename.split('_')
                    if 'PIRADS' in filename_list and filename_list[len(filename_list)-1]=='bt.voi' and not filename_list[1]=='p':
                        if not filename_list[4]=='PIRADS':
                            problems+=[root+'_'+filename+'_'+' location issue']
        print(sorted(problems))

#################################################################################################################

class DatasetSummary:

    def __init__(self):
        self.basePATH=r'C:\Users\sanfordt\Desktop\PIRADS_dataset\databases'
        self.database='prostateX_new'
        self.workingdb='revision_analysis_2'
        self.devfolder='model_dev_indvPIRADS'

    def lesions_summary(self):
        '''in file called 'tumors' evaluates the dataset characteristics'''
        tumors=os.listdir(os.path.join(self.basePATH,self.database,'tumors')) #total number of tumors
        studies=pd.Series(file.split('_')[0] for file in os.listdir(os.path.join(self.basePATH,self.database,'lesions_by_study'))) #total number MRIs
        duplicated=studies[studies.duplicated()] #check to make sure there are no duplicates
        print("The following studies are duplicated {}".format(duplicated))
        num_unique_patients=len(list(set([file.split('_')[0] for file in os.listdir(os.path.join(self.basePATH,self.database,'tumors'))])))

        side=[]; location=[]; zone=[]; PIRADS=[];
        for tumor in tumors:
            side += [tumor.split("_")[3]]
            location += [tumor.split("_")[4]]
            zone += [tumor.split("_")[5]]
            PIRADS += [tumor.split("_")[7]]


        print("number tumors is is {}".format(len(tumors)))
        print("number of studies is {}".format(studies.count()))
        print("number of patients is {}".format(num_unique_patients))
        print("side of right sided tumors is  {}".format(Counter(side)))
        print('location of patients is {}'.format(Counter(location)))
        print('zone of patients is {}'.format(Counter(zone)))
        print('PIRADS score is {}'.format(Counter(PIRADS)))

    def slice

    def demographic_data(self,db='consecutive'):
        '''
        obtaining demographics for each dataset individually.
        '''

        patients=os.listdir(os.path.join(self.basePATH,self.database,db)) #lesions with PIRADS 2 or greater
        tumors_included=os.listdir(os.path.join(self.basePATH,self.workingdb,'lesions_by_study'))  #change this part if you update your database, this is all patients included in the database
        excluded=set(patients).difference(set(tumors_included))  #patients that have PIRADS 1
        patients_overlap=list(set(patients).intersection(set(tumors_included))) #patient that have PIRADS 2 or greater
        mrns=pd.Series([patient.split('_')[0] for patient in patients_overlap])  #get all mrns in your database
        duplicated=mrns[mrns.duplicated()] #sanity check to make sure no repeat patients.

        #loop over patients, read in dicom files and extract demographic data
        age_list=[]; weight_list=[]; total=0
        for patient in patients_overlap:
            t2_path=os.path.join(self.basePATH,self.database,db,patient,'dicoms','t2')
            dcm=pydicom.dcmread(os.path.join(t2_path,os.listdir(t2_path)[0]))
            weight=dcm[0x00101030].value
            DOB = dcm[0x00100030].value
            DOB_datetime = date(year=int(DOB[0:4]), month=int(DOB[4:6]), day=int(DOB[6:8]))
            age = calculate_age(DOB_datetime)  #helper function
            age_list+=[age]
            weight_list+=[weight]
            total+=1

        weight_list=[val for val in weight_list if val>18]
        median_age=statistics.median(age_list)
        min_age=min(age_list); max_age=max(age_list)
        median_weight=statistics.median(weight_list)
        min_weight=min(weight_list); max_weight=max(weight_list)

        print("for database {}".format(db))
        print("total of {} studies".format(len(patients)))
        print("total of {} studies excluded".format(len(excluded)))
        print("The following studies are duplicated {}".format(duplicated))
        print("total of {} studies with >PIRADS 2".format(len(patients_overlap)))
        print("median age is {} with min age of {} and max age of {}".format(median_age,min_age,max_age))
        print("median weight is {} with min weight of {} and max weight of {}".format(median_weight,min_weight,max_weight))

    def ER_coil(self, db='consecutive'):
        '''this function returns a count of the highB value headers.
        At our instiution, b value of 2000 used only with endorectal coil.  All others are with b1500'''
        out_list = []
        patients = os.listdir(os.path.join(self.basePATH, self.database, db))
        tumors_included = os.listdir(os.path.join(self.basePATH, self.workingdb,'lesions_by_study'))  # change this part if you update your database
        patients_overlap = list(set(patients).intersection(set(tumors_included)))
        for patient in patients_overlap:
            highbs = os.listdir(os.path.join(self.basePATH,self.database,db,patient,'dicoms', 'highb','raw'))
            ds = pydicom.dcmread(os.path.join(self.basePATH,self.database,db,patient,'dicoms', 'highb','raw', highbs[0]))
            out_list += [ds[0x08, 0x103e].value]
        print(Counter(out_list))

    def calc_volumes(self,filetype='wp'):
        '''
        calculate the volumes for a segmented structure based on .nifti mask volume
        :return:
        '''

        base_path=os.path.join(self.basePATH,self.database)
        outDF=pd.DataFrame()
        for patient in os.listdir(base_path):
            # calculate voxel size size
            first_t2=os.listdir(os.path.join(base_path, patient, 'dicoms', 't2'))[0]
            ds=pydicom.dcmread(os.path.join(base_path,patient,'dicoms','t2',first_t2))
            xy_size=ds[0x28,0x30].value; z_size=ds[0x18,0x88].value
            volume_voxel=xy_size[0]*xy_size[1]*z_size

            #search the files for the type you are interested in (i.e. wp=whole prostate)
            filelist = []
            for file in os.listdir(os.path.join(base_path, patient,'nifti','mask')):
                if len(file.split('_'))<5:
                    if filetype == 'wp': pat = re.compile('([Ww][Pp]){1}')
                    elif filetype == 'tz': pat = re.compile('([Tt][Zz]){1}')
                    if re.search(pat,file) !=None: filelist+=[file]

            #select annotations for expert
            for i in range(len(filelist)):
                name=filelist[i]
                name_noend=name.split('.nii')[0]
                if name_noend.split('_')[-1]=='bt':filename=name
                elif name_noend.split('_')[-1]=='mm':filename=name
                elif name_noend.split('_')[-1]=='ts':filename=name
                elif name_noend.split('_')[-1] == 'pseg': filename = name
                elif name_noend.split('_')[-1] == 'dk':filename = name

            #calculate volume
            nifti_path=os.path.join(base_path,patient,'nifti','mask',name)
            volume=calculate_volume(nifti_path, patient, volume_voxel)  #helper function listed below
            print("volume of {} for patient {} is: {}".format(filetype,patient.split('_')[0],volume))
            series = pd.DataFrame([patient, volume]).transpose()
            outDF = pd.concat([outDF, series], axis=0)

        #save all volumes to file
        outDF.to_csv(os.path.join(self.basePATH,self.database+'_volumes_of_'+filetype+'.csv'))


###########################################################################
####################### Helper Functions #################################
###########################################################################


def get_slices(path='',dir='', ext='.jpg'):
    '''recursively looks in folder for files with specific file extension'''
    num = 0
    list_filenames = []
    print("looking for all slices in directory {}".format(dir))
    for root, dirnames, filenames in os.walk(os.path.join(path,dir)):
        for filename in filenames:
            if filename.endswith(ext):
                num += 1
                list_filenames += [filename]
    print("total of {} files in the directory {} ".format(len(list_filenames), dir))
    return (list_filenames)


def calculate_volume(path,patient, volume_voxel):
    # calculate volume
    volume = nib.load(path)
    vol_array = volume.get_fdata()
    volume = (int(vol_array.sum()) * volume_voxel) / 1000
    volume = round(volume, 2)
    return volume

def calculate_age(born):
    today = date.today()
    return today.year - born.year - ((today.month, today.day) < (born.month, born.day))

if __name__=='__main__':
    c=DatasetInterrogate()
    c.check_pt_spillage()