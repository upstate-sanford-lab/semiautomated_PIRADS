import os
import sys
import datetime

import fastai.vision as faiv
import numpy as np
import matplotlib.pyplot as plt
import torch
from torch import nn
import PIL
from PIL import Image
from resizeimage import resizeimage

from fastai.vision import *
from fastai.callbacks import *
from fastai.callbacks.tracker import *
from fastai.basic_train import *

#author @t_sanf

class ModelApply:


    def __init__(self):
        self.imagedir = '/home/tom/Desktop/revision_analysis2/model_dev_indvPIRADS'
        self.outdir = '/home/tom/Desktop/revision_analysis2/output'
        self.testPath = os.path.join(self.imagedir, 'val_test_pt') #need to make new directory with all validation and test samples
        self.clin_val ='/home/tom/Desktop/revision_analysis2/annotation'
        self.save_dir= os.path.join(self.outdir,'eval_model')
        self.device=0


    def convert_model_to_export(self, model_name):
        initial_filename = os.path.join(self.outdir, 'exported_models', model_name)
        final_filename = os.path.join(self.outdir, 'exported_models', 'export.pkl')
        shutil.copy2(initial_filename, final_filename)

    def assign_val_test(self):
        '''assigns 'test' or 'val' based on file structure to file called 'tumor_level_annoation' and re-assigns to name 'clinical validation' '''

        data=pd.read_csv(os.path.join(self.clin_val,'tumor_level_annotation.csv'))
        test_data=os.listdir(os.path.join(self.imagedir,'test'))
        val_data=os.listdir(os.path.join(self.imagedir,'val_pt'))

        include=0
        for name in data['tumor_name'].tolist():
            if name in test_data:
                include+=1
                data.loc[data['tumor_name']==name,'val_test']='test'
            if name in val_data:
                include+=1
                data.loc[data['tumor_name']==name,'val_test']='val'


        data_u=data.dropna(subset=['val_test'])
        data_u.to_csv(os.path.join(self.clin_val,'clinical_validation.csv'))

    def apply_test_vote(self):
        torch.cuda.set_device(self.device)

        # set up the output directory
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        if not os.path.isdir(os.path.join(self.outdir, 'val_results')):
            os.mkdir(os.path.join(self.outdir, 'val_results'))

        #import model
        model_path = os.path.join(self.outdir, 'exported_models')
        learn = load_learner(model_path)
        print("learner loaded")
        test_path = self.testPath

        #import clinical data
        tumor_val_db=pd.read_csv(os.path.join(self.clin_val,'clinical_validation.csv'))
        df_out = pd.DataFrame()

        print('making predictions')
        for tumor in os.listdir(os.path.join(test_path)):
            sum_pred = np.zeros(4)
            square_pred = np.zeros(4)
            img_num = 0

            for image in sorted(os.listdir(os.path.join(test_path, tumor))):
                img = open_image(os.path.join(test_path, tumor, image),resize=False)
                pred_class, pred_idx, outputs = learn.predict(img)
                sum_pred += outputs.numpy()
                square_pred += (outputs.numpy()) ** 2
                img_num += 1

            # metrics
            average = sum_pred / img_num
            ave_pred_class = np.argmax(average)

            # make prediction change this part to change the method of classification
            pred_PIRADS = str(ave_pred_class + 2)
            gt_PIRADS = image.split('_')[8]

            rowIndex = tumor_val_db.index[tumor_val_db.loc[:,'tumor_name']==tumor]
            tumor_val_db.loc[rowIndex, 'pred_PIRADS'] = pred_PIRADS
            tumor_val_db.loc[rowIndex, 'gt_PIRADS'] = gt_PIRADS

        tumor_val_db.to_csv(os.path.join(self.save_dir,'tumor_val_db.csv'))


    def PIRADS_mapping(self,preprocess=False,col='gt_PIRADS'):
        '''figure out how well PIRADS maps to outcome'''

        if preprocess==True:
            self.apply_test_vote()

        db=pd.read_csv(os.path.join(self.save_dir,'tumor_val_db.csv'))
        db_bx=db.loc[db.loc[:,'biopsy']==1]
        total_tumors=db_bx.shape[0]
        print("total tumors {}".format(total_tumors))

        for score in [2,3,4,5]:
            print(score)
            CA = 0
            CS_CA = 0
            total_t_P = 0
            agreement=0
            within_1=0
            df_P = pd.DataFrame(columns=db_bx.columns)

            for tumor in db_bx.loc[:,'tumor_name']:
                t_series=db_bx.loc[db_bx.loc[:,'tumor_name']==tumor,:]
                if int(t_series[col])==int(score):
                    df_P=df_P.append(t_series)

            df_P.to_csv(os.path.join(self.save_dir,'tumor_val_'+str(score)+'.csv'))

            #determine mapping to clinical significance
            for tumor in df_P.loc[:,'tumor_name']:

                #find gleason
                total_t_P += 1
                p_series = df_P.loc[df_P.loc[:, 'tumor_name'] == tumor, :]
                primary_G=int(p_series.loc[:,'primary_gleason'])
                sec_G=int(p_series.loc[:,'seconary_gleason'])
                total_G=primary_G+sec_G
                if total_G>0:
                    CA+=1
                if total_G>=7:
                    CS_CA+=1

                #find PIRADS
                gt_PIRADS=int(p_series.loc[:,'gt_PIRADS'])
                pred_PIRADS=int(p_series.loc[:,'pred_PIRADS'])
                if gt_PIRADS==pred_PIRADS:
                    agreement+=1
                if abs(gt_PIRADS-pred_PIRADS)<2:
                    within_1+=1


            print("for PIRADS {}, of a total of {} patients:".format(score,total_t_P))
            print("cancer detection rate is {}".format(round(CA/total_t_P*100,2)))
            print("clinically significant detection rate is {}".format(round(CS_CA/total_t_P*100,2)))
            print("agreement in {} of {}".format(agreement,total_t_P))
            print("agreement rate is {}".format(round(agreement/total_t_P*100,2)))
            print("within 1 rate in {} of {}".format(within_1,total_t_P))
            print("within 1 rate is {}".format(round(within_1/total_t_P*100,2)))
            print("------------")
        print("total tumors {}".format(total_tumors))



    def PIRADS_agreement(self,preprocess=True,db_train_test=True,db_val='test',col='gt_PIRADS'):
        '''figure out how well PIRADS maps to outcome'''

        if preprocess==True:
            self.apply_test_vote()

        db=pd.read_csv(os.path.join(self.save_dir,'tumor_val_db.csv'))


        if db_train_test==True:
            if db_val=='val':
                db=db.loc[db.loc[:,'val_test']=='val',:]
            if db_val=='test':
                db = db.loc[db.loc[:, 'val_test'] == 'test', :]

        total_tumors = db.shape[0]

        total_agreement=0
        total_within_1=0
        for score in [2,3,4,5]:
            print(score)
            agreement=0
            within_1=0
            total_t_P =0
            df_P = pd.DataFrame(columns=db.columns)


            for tumor in db.loc[:,'tumor_name']:
                t_series=db.loc[db.loc[:,'tumor_name']==tumor,:]
                if int(t_series[col])==int(score):
                    df_P=df_P.append(t_series)

            #determine mapping to clinical significance
            for tumor in df_P.loc[:,'tumor_name']:

                #find PIRADS
                total_t_P += 1
                p_series = db.loc[db.loc[:, 'tumor_name'] == tumor, :]

                #find PIRADS
                gt_PIRADS=int(p_series.loc[:,'gt_PIRADS'])
                pred_PIRADS=int(p_series.loc[:,'pred_PIRADS'])
                if gt_PIRADS==pred_PIRADS:
                    agreement+=1
                if abs(gt_PIRADS-pred_PIRADS)<2:
                    within_1+=1
            total_agreement+=agreement
            total_within_1+=within_1

            print("agreement in {} of {}".format(agreement,total_t_P))
            print("agreement rate is {}".format(round(agreement/total_t_P*100,2)))
            print("within 1 rate in {} of {}".format(within_1,total_t_P))
            print("within 1 rate is {}".format(round(within_1/total_t_P*100,2)))
            print("------------")
        print("total tumors {}".format(total_tumors))
        print("total agreement {}".format(total_agreement))
        print("total agreement rate is {}".format(round(total_agreement/total_tumors*100,2)))
        print("within 1 agreement {}".format(total_within_1))
        print("total within 1 agreement is {}".format(round(total_within_1/total_tumors*100,2)))



if __name__ == '__main__':
    c = ModelApply()
    c.convert_model_to_export(model_name='model_name.pkl')
    #c.PIRADS_agreement()
    #c.PIRADS_mapping(preprocess=False, col='pred_PIRADS')
