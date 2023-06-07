__cmddoc__ = """

drugpp.py - Reads an SDFile or SMIFIle molecules with Synthetic Accessibility score using RDKIT tools
#Contact: haiming_cai@hotmail.com - 2022 - CHINA-VMI

ps -ef | grep 'Google Chrome Helper' | awk '{print $2}' | xargs kill -9
conda activate drugfilter_py39_env
python proxy.drug.admet.fix7.py refsmile.smi Cal_ademt.csv

""" 
import os
import io
import csv
import sys
import requests
from urllib.parse import urljoin
from bs4 import BeautifulSoup
from fake_useragent import UserAgent

import pandas as pd
import numpy as np
import time

from rdkit import Chem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import requests
from typing import Tuple
import tkinter
import random
import asyncio
#import traceback
import pyppeteer
from pyppeteer import launch
from pyppeteer.launcher import Launcher
' '.join(Launcher().cmd)

import warnings
warnings.filterwarnings('ignore')

def get_proxy():
    return requests.get("http://127.0.0.1:5010/get/").json()

def delete_proxy(proxy):
    requests.get("http://127.0.0.1:5010/delete/?proxy={}".format(proxy))

def sdf2smi(sdflist, key1):
     mylogger = RDLogger.logger()
     mylogger.setLevel(val=RDLogger.ERROR)

     #suppl = []
     #for mol in Chem.SDMolSupplier(sdflist)
     #    if mol:
     #       try:
     #          suppl.append(mol)
     #       except:
     #          pass 

     suppl = [ mol for mol in Chem.SDMolSupplier(sdflist) ]

     ID_list = []
     smiles_list = []
     for mol in suppl:
         if mol:
             name = mol.GetProp(key1)
             ID_list.append(name)
             try:
                smi = Chem.MolToSmiles(mol,isomericSmiles=False)
                smiles_list.append(smi)
             except:
                pass
     mol_datasets = pd.DataFrame({'smiles':smiles_list,'ID':ID_list})
     #mol_datasets.to_csv(os.path.join('./',"TMDB.smi"), sep='\t', index=False, header=False, encoding='utf-8')
     mol_datasets = mol_datasets[['smiles','ID']]
     return mol_datasets

def smi2list(smifile, key="_Name"):
     #mylogger = RDLogger.logger()
     #mylogger.setLevel(val=RDLogger.ERROR)

     #print(smifile)
     suppl = Chem.SmilesMolSupplier(smifile, delimiter='\t', titleLine=False, nameColumn=1)
     ID_list = []
     smiles_list = []
     #print(suppl)
     for mol in suppl:
        if mol:
           try:
              smi = Chem.MolToSmiles(mol,isomericSmiles=False)
              #print(smi)
              smiles_list.append(smi)
              name = mol.GetProp(key)
              #print(name)
              ID_list.append(name) 
           except:
              pass

     mol_datasets = pd.DataFrame({'smiles':smiles_list,'ID':ID_list})
     mol_datasets = mol_datasets[['smiles','ID']]
     return mol_datasets

def canonic_smiles(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    return Chem.MolToSmiles(mol)

def swissadme(smidf, outfile = 'Cal_adme.csv'):
    root_url = "http://www.swissadme.ch/"
    url = "http://www.swissadme.ch/index.php"
    totaldf = pd.DataFrame()
    ua = UserAgent()
    headers = {
        'User-Agent': ua.random
    }
    proxy = get_proxy().get("proxy")
    file2 = ""
    for index,value in smidf.iterrows():
        name = value['ID']
        smiles = value['smiles']
        print(name, smiles)
        idresult = {"ID":name,"Origin smiles":smiles}
        iddf = pd.DataFrame(idresult,index=[0])
        csmiles =  canonic_smiles(smiles) 
        data = {"smiles": csmiles}
        for i in range(1, 101):
            try:
                print("The request retry time is %d" % i)
                requests.adapters.DEFAULT_RETRIES = 2
                sess = requests.Session()
                sess.keep_alive = False
                r1 = sess.post(url, data = data, headers=headers, proxies={"http": "http://{}".format(proxy)}, verify=False, timeout=(100,200))
                if r1.status_code == 200:
                   try:
                      soup = BeautifulSoup(r1.text, features = 'html.parser')
                      for a in soup.select('a[href^="results"]'):
                         file_url = a['href']
                      file = urljoin(root_url, file_url)
                   except Exception:
                      raise Exception
            
                   if file == file2:
                      secends = random.randint(5, 40)
                      time.sleep(secends)
                      raise Exception
                   else:
                      try:
                         #requests.adapters.DEFAULT_RETRIES = 2
                         #sess = requests.Session()
                         #sess.keep_alive = False
                         r2 = sess.get(file, proxies={"http": "http://{}".format(proxy)}, headers=headers, verify=False, timeout=(100, 200))
                      except Exception:
                         raise Exception

                      if r2.status_code == 200:
                         file_data = r2.content
                      else:
                         raise Exception
                else: 
                   raise Exception

                sess.close()
                secends = random.randint(5, 40)
                time.sleep(secends)
                break
            except Exception:
                sess.close()
                secends = random.randint(5, 40)
                time.sleep(secends)
                headers = {
                    'User-Agent': ua.random
                }
                proxy = get_proxy().get("proxy")
                continue

        file2 = file    

        resultdf = pd.read_csv(io.StringIO(file_data.decode('utf-8')))
        #print(resultdf)
        resultdf = pd.concat([iddf,resultdf.reset_index(drop=True)],axis=1)
        totaldf = pd.concat([totaldf,resultdf],axis=0,ignore_index=False)
    totaldf.drop(['Molecule'], axis=1,inplace=True)
    #totaldf = totaldf.rename(columns={"Canonical smiles": "SMILES"})
    totaldf = totaldf.rename(columns = lambda x: x.replace("#",""))
    #totaldf.to_csv(outfile,index=False)
    return totaldf

class Logger:
    def __init__(self, filepath) -> None:
        self.filepath = filepath

    def write_log(self, txt):
        with open(self.filepath, 'a') as f:
            f.write(f"{txt}\n")
            print(txt)

async def antiAntiCrawler(page):
    await page.evaluate('''() =>{ Object.defineProperties(navigator,{ webdriver:{ get: () => undefined } }) }''')
    await page.evaluate('''() =>{ window.navigator.chrome = { runtime: {},  }; }''')
    await page.evaluate('''() =>{ Object.defineProperty(navigator, 'languages', { get: () => ['en-US', 'en'] }); }''')
    await page.evaluate('''() =>{ Object.defineProperty(navigator, 'plugins', { get: () => [1, 2, 3, 4, 5,6], }); }''')	

async def navigate_pkcsm_site(page, drug_name, canonical_smile, logger): 
    #print('two')

    logger.write_log("Navigating pkcsm site...")

    pkcsm_url = 'https://biosig.lab.uq.edu.au/pkcsm/prediction'
    await page.goto(pkcsm_url, {"timeout": 100000, "waitUntil": "domcontentloaded"})

    pkcsm_js_code = f"""
        smile_input = document.querySelector(
            `body > div.container > div.row > div.span7.offset1 
            > form > div:nth-child(2) > div:nth-child(3) > div > input`
            );        
        smile_input.value = "{canonical_smile}"
    """
    await page.evaluate(pkcsm_js_code)

    navpromise = asyncio.ensure_future(page.waitForNavigation({"timeout": 60000, "waitUntil": "networkidle0"})) 

    await page.click("body > div.container > div.row > div.span7.offset1 >" + 
            "form > div:nth-child(5) > div > div > div > div:nth-child(2) > button")

    await navpromise

    logger.write_log("Waiting for pkcsm calculation...")
    pkcsm_result =  await wait_till_computation(page, drug_name)

    logger.write_log("Scrapped pkcsm site...")
    
    return pkcsm_result

async def wait_till_computation(page, drug_name) -> Tuple:
    #print('one')
    content = await page.content()
    soup = BeautifulSoup(content, "html.parser")
    nextRequest = asyncio.ensure_future(page.waitForNavigation({"timeout": 6000, "waitUntil": "networkidle0"})) 
    
    is_all_loaded = True

    try:
        await nextRequest 
    except pyppeteer.errors.TimeoutError: # we will consider that all of our data is loaded
        pass

    data_table = soup.select_one("body > div.container > div.row.fluid > div.span8 > div.well > table > tbody")
    all_props = data_table.find_all("tr")

    # check if calculation is still running
    for prop in all_props:
        table_data = prop.find_all("td")
        if("Running" in table_data[2].text):
            is_all_loaded = False
    
    if not is_all_loaded: 
        result = await wait_till_computation(page, drug_name) 
        return result
    
    else: 
        return scrape_results(soup, drug_name)

def scrape_results(soup, drug_name) -> Tuple:
    #print('three')

    data_table = soup.select_one("body > div.container > div.row.fluid > div.span8 > div.well > table > tbody")
    all_props = data_table.find_all("tr")

    #water_solubility = all_props[0].find_all("td")[2].text.strip()
    #intestinal_absorption = all_props[2].find_all("td")[2].text.strip()
    #bbb_permeability = all_props[9].find_all("td")[2].text.strip()
    #cns_permeability = all_props[10].find_all("td")[2].text.strip()
    #cyp3a4_substrate = all_props[12].find_all("td")[2].text.strip()
    #cyp3a4_inhibitor = all_props[17].find_all("td")[2].text.strip()
    #total_clearance = all_props[18].find_all("td")[2].text.strip()
    #renal_oct2_substrate = all_props[19].find_all("td")[2].text.strip()
    ames_toxicity = all_props[20].find_all("td")[2].text.strip()
    max_tolerated_dose = all_props[21].find_all("td")[2].text.strip()
    hERG_I_inhibitor = all_props[22].find_all("td")[2].text.strip()
    hERG_II_inhibitor = all_props[23].find_all("td")[2].text.strip()    
    oral_rat_acute_toxicity = all_props[24].find_all("td")[2].text.strip()
    hepatotoxicity = all_props[26].find_all("td")[2].text.strip()
    skin_sensitisation = all_props[27].find_all("td")[2].text.strip()
    minnow_toxicity = all_props[29].find_all("td")[2].text.strip()  

    #absorption_col = (water_solubility, intestinal_absorption)
    #distribution_col = (bbb_permeability, cns_permeability)
    #metabolism_col = (cyp3a4_substrate, cyp3a4_inhibitor)
    #excretion_col = (total_clearance, renal_oct2_substrate)
    #toxicity_col = (max_tolerated_dose, oral_rat_acute_toxicity)

    #total_row = (drug_name,) + absorption_col + distribution_col + metabolism_col + excretion_col + toxicity_col
    total_row = [ames_toxicity, max_tolerated_dose, hERG_I_inhibitor, hERG_II_inhibitor, oral_rat_acute_toxicity, hepatotoxicity, skin_sensitisation, minnow_toxicity]
    return total_row

async def pkcsm(smidf, outfile = 'Cal_toxity.csv'):
    head = ["Ames toxicity", \
            "Max tolerated dose", \
            "hERG I inhibitor", \
            "hERG II inhibitor", \
            "Oral rat acute toxicity", \
            "Hepatotoxicity", \
            "Skin sensitisation", \
            "Minnow toxicity"]

    totaldf = pd.DataFrame()
    logger = Logger('./log.log')

    ua = UserAgent()
    #agent = ua.random
    agent = 'Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.67 Safari/537.36'
    #proxy = "http://"+get_proxy().get("proxy")
    proxy = get_proxy().get("proxy")
    #print(proxy)
    for index,value in smidf.iterrows():
        name = value['ID']
        smiles = value['smiles']
        print(name, smiles)
        idresult = {"ID":name,"Origin smiles":smiles}
        iddf = pd.DataFrame(idresult,index=[0])
        csmiles =  canonic_smiles(smiles) 

        for i in range(1, 101):
            try:
                print("The request retry time is %d" % i)
                browser = await launch(headless=True, args=['--no-sandbox', '--disable-dev-shm-usage', '--disable-gpu', '--no-zygote', '--ignore-certificate-errors', '--enable-features=NetworkService', "--proxy-server={}".format(proxy)])
                page = await browser.newPage()
                #time.sleep(10)
                await page.setUserAgent('Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.67 Safari/537.36')
                await page.setViewport({'width':1920, 'height':1080})
                pkcsm_result = await navigate_pkcsm_site(page, name, csmiles, logger)
            #except (pyppeteer.errors.PageError, pyppeteer.errors.TimeoutError) as e:
            except Exception as e:
                #print(e)
                await browser.close()
                time.sleep(2)
                #delete_proxy(proxy)
                code = 403
                while code != 200:
                    proxy = get_proxy().get("proxy")
                    pkcsm_url = 'https://biosig.lab.uq.edu.au/pkcsm/prediction'
                    #print(code)
                    try:
                        browser = await launch(headless=True, args=['--no-sandbox', '--disable-dev-shm-usage', '--disable-gpu', '--no-zygote', '--ignore-certificate-errors', '--enable-features=NetworkService', "--proxy-server={}".format(proxy)])
                        page = await browser.newPage()
                        #await page.setUserAgent(agent)
                        await page.setUserAgent('Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.67 Safari/537.36')
                        await page.setViewport({'width':1920, 'height':1080})
                        res = await page.goto(pkcsm_url, options={"timeout": 30000})
                        code = res.status
                        #print(code)
                        await asyncio.sleep(2)
                        print("Get a success proxy: {}".format(proxy))
                    except Exception:
                        print("The proxy {} is failed.".format(proxy))
                        await browser.close()
                        #delete_proxy(proxy)
                        code = 403
                #print(proxy)
                await browser.close()
                time.sleep(10)
                continue
            else:
                await browser.close()
                time.sleep(30)
                break

        #print(pkcsm_result)
        #print(len(pkcsm_result))
        #print(len(head))
        pkcsm_details = pd.DataFrame([pkcsm_result], columns = head)
        #print(pkcsm_details)
        resultdf = pd.concat([iddf, pkcsm_details.reset_index(drop=True)], axis=1)
        #print(resultdf)
        totaldf = pd.concat([totaldf, resultdf], axis=0, ignore_index=False)
        #print(totaldf)
        #await browser.close()
    #totaldf.to_csv(outfile,index=False)
    #print(totaldf)
    return totaldf
 
if __name__ == '__main__':
    import sys
    smifile = sys.argv[1]
    outfile = sys.argv[2]

    start = time.time()
    print("start time:", time.asctime(time.localtime(start)))

    if os.path.getsize(smifile) != 0:
       dataset1 = smi2list(smifile)
       dataset2 = swissadme(dataset1)
       dataset2 = dataset2[dataset2.columns.drop(list(dataset2.filter(regex='Unnamed')))]
       dataset3 = asyncio.get_event_loop().run_until_complete(pkcsm(dataset1))
       dataset4 = pd.merge(dataset2, dataset3, on=['ID', 'Origin smiles'], how='outer')
       dataset4.to_csv(outfile,index=False)
    else:
       print("The input smiles file may be null or wrong, please check it again.")

    end = time.time()
    print("The %s job finished without error!\nTime is (s):" % smifile.split("/")[-1].split("_filter")[0], end - start)
