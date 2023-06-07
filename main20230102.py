import sys
import os
import io
import pandas as pd
from rdkit import rdBase
from rdkit import Chem
import argparse
import logging
import time

logger = logging.getLogger()
logger.level = logging.INFO
stream_handler = logging.StreamHandler(sys.stderr)
logger.addHandler(stream_handler)
os.environ['NUMEXPR_MAX_THREADS'] = '16'

def admetfilter(admetdf, outfile = 'admet.filter.smi'):
    filterdf = admetdf.query("(`BBB permeant` == 'No') & \
                              (`GI absorption` == 'High') & \
                              (`Pgp substrate` == 'No') & \
                              (`CYP1A2 inhibitor` == 'No') & \
                              (`CYP2C19 inhibitor` == 'No') & \
                              (`CYP2C9 inhibitor` == 'No') & \
                              (`CYP2D6 inhibitor` == 'No') & \
                              (`CYP3A4 inhibitor` == 'No') & \
                              (`Lipinski violations` == 0 or `Lipinski violations` == 1) & \
                              (`Veber violations` == 0) & \
                              (`Ames toxicity` == 'Yes') & \
                              (`hERG I inhibitor` == 'No') & \
                              (`hERG II inhibitor` == 'No') & \
                              (`Hepatotoxicity` == 'No') & \
                              (`Skin sensitisation` == 'No') & \
                              (`PAINS alerts` == 0)")
    filtersmidf = filterdf[['Origin smiles', 'ID']]
    filtersmidf.to_csv(outfile, sep='\t', index=False, header=False, encoding='utf-8')

def clusterdf2admet():
    os.system('cp ./src/admet2.fix3.sh ./admet2.sh')
    for file in os.listdir('./'):
        if file.endswith('.smi'):
           targetname = file.split(".")[0]
           filepath = "./" + file
           for i in range(1, 4):
               try:
                  os.system('./admet2.sh %s' % filepath) 
                  dataset = pd.read_csv("Cal_ademt.csv")
                  admetfilter(dataset)
                  os.system('mv Cal_ademt.csv ./%s_Cal_ademt.csv' % targetname)
                  os.system('mv admet.filter.smi ./%s_admet.filter.smi' % targetname)  
                  os.system('rm -rf log.log') 
                  break
               except Exception as e:
                  os.system('./src/update.proxy.sh')           
                  os.system('rm -rf liglist* Cal_ademt.csv log.log')             
    os.system('rm -rf ./admet2.sh')
    os.system('mkdir -p ./ok')

def main_loop():
    if os.path.exists('./ok') == False:
       logging.info('The admet job had been upload!')
       clusterdf2admet()
       logging.info('The admet job had been finished!')
    else:
       logging.info('The admet job had been finished!')

if __name__ == '__main__':
    main_loop()  