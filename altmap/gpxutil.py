# standard libs
import urllib.parse
import re
import os
import json
import sys
import gc
import pkgutil
import pathlib

# external libs
import requests
from tqdm import tqdm
import numpy as np
import gpxpy
import ndjson

def getdatapath():
    package = pkgutil.get_loader("altmap")
    datapath = str(pathlib.Path(package.get_filename()).parent)+"/data/"

    return datapath


def getGpxurl(keyword,placename=None,starttime=None,endtime=None,area=None,genre=None,maxgpxnum=30):

    query="pnum=1&"

    if placename != None:
        qplacename = "mname="+urllib.parse.quote(keyword,encoding='euc-jp')+"&"
        query += qplacename
    else:
        query += "mname=&"

    qkeyword = "place="+urllib.parse.quote(keyword,encoding='euc-jp')+"&"
    query += qkeyword

    if starttime != None:
        qstarttime = "start_y="+starttime.split("-")[0]+"&start_m="+starttime.split("-")[1]+"&start_d="+starttime.split("-")[2]+"&"
        query += qstarttime
    else:
        query += "start_y=&start_m=&start_d=&"

    if endtime != None:
        qendtime   = "end_y="+endtime.split("-")[0]  +"&end_m="+endtime.split("-")[1]  +"&end_d="+endtime.split("-")[2]+"&"
        query += qendtime
    else:
        query += "end_y=&end_m=&end_d=&"

    if area !=None:
        qarea = "area="+str(int(area))+"&"
        query += qarea
    else:
        query += "area=&"

    if genre !=None:
        qgenre = "genre="+str(int(genre))+"&"
        query += qgenre
    else:
        query += "genre=&"
    
    qplacename = "mname="+urllib.parse.quote(keyword,encoding='euc-jp')+"&"
    query += qplacename

    qkeyword = "place="+urllib.parse.quote(keyword,encoding='euc-jp')+"&"
    query += qkeyword

    query += "istrack=1&uname=&order=start&request=1&relpts=&submitbtn=%B8%A1%BA%F7"

    url = "https://www.yamareco.com/modules/yamareco/search_record.php?"+query
    res = requests.get(url)
    tmp = res.content
    tmp.decode('euc-jp')
    contents = res.content.decode('euc-jp')
    findres = np.unique(re.findall("https://www\.yamareco\.com/modules/yamareco/detail-[0-9]*[0-9]\.html",contents))
    
    gpxurlList = []
    for url in findres:
        gpxurlList.append( re.sub("html", "gpx", re.sub("detail", "track", url) ) )
    
    return gpxurlList[0:min(len(gpxurlList),maxgpxnum)]

def loadGpxfromurl(gpxurl):
    
    urldata = requests.get(gpxurl).content
    gpx = gpxpy.parse(urldata)
    
    return gpx


def loadDatafromndjson(path):

    with open(path,"r") as f:
        input_data = ndjson.load(f)
    
    return input_data

def getArealist():
    return loadDatafromndjson(getdatapath()+"area.ndjson")

def getGenrelist():
    return loadDatafromndjson(getdatapath()+"genre.ndjson")

def getTypeidlist():
    return loadDatafromndjson(getdatapath()+"typeid.ndjson")

def searchPointdata(
    coordinate_search=True, lat_max=0, lat_min=0, lon_max=0, lon_min=0,
    typeid_search=False, typeid="1",   
    ):

    temp = loadDatafromndjson(getdatapath()+"pointdata.ndjson")

    reslist = [d for d in temp
    if (
        (~coordinate_search or ( (lat_min < float(d['lat'])) and (lat_max > float(d['lat'])) and (lon_min < float(d['lon'])) and (lon_max > float(d['lon'])) )) and
        (~typeid_search or ( (d.get('types') !=None) and (typeid in [temptype.get('type_id') for temptype in d.get('types')]) ))
    )]

    return reslist

def savePointdatafromAPI(savepath=None,pagemin=1,pagemax=800,overwrite=False):

    if savepath==None:
        _savepath = getdatapath()+"pointdata.ndjson"
    
    else:
        _savepath = savepath
    
    if overwrite:
        os.system("rm -rf "+_savepath)

    for i in tqdm(range(pagemax)):
        npage=i+pagemin
        typeid=0
        url = "https://api.yamareco.com/api/v1/searchPoi"
        params = {"page":npage,"type_id":typeid}
        res = requests.post(url,data=params)
        json_str = res.content.decode("euc-jp")
        json_dict = json.loads(json_str)

        flag = False
        try:
            temp = json_dict['poilist'][0]['ptid']
            flag = True
        except:
            break

        if flag:
            with open(_savepath,"a") as f:
                writer = ndjson.writer(f)
                [writer.writerow(comp) for comp in json_dict['poilist']]

def updatePointdatafromAPI(savepath=None,pagemax=900):

    if savepath==None:
        _savepath = getdatapath()+"pointdata.ndjson"
    
    else:
        _savepath = savepath
    
    if not pathlib.Path(_savepath).exists():
        print('Error: Pointdata file "'+_savepath+'" dose not exist', file=sys.stderr)
        sys.exit(1)
    
    with open(_savepath,"r") as f:
        temp = ndjson.load(f)
        lastptid = int(temp[-1]["ptid"])
        del temp
        gc.collect
    
    for i in tqdm(range(pagemax-780)):
        npage=i+780
        typeid=0
        url = "https://api.yamareco.com/api/v1/searchPoi"
        params = {"page":npage,"type_id":typeid}
        res = requests.post(url,data=params)
        json_str = res.content.decode("euc-jp")
        json_dict = json.loads(json_str)

        flag = False
        try:
            temp = json_dict['poilist'][0]['ptid']
            flag = True
        except:
            break
        
        if flag:
            with open(_savepath,"a") as f:
                writer = ndjson.writer(f)
                [writer.writerow(comp) for comp in json_dict['poilist'] if int(comp['ptid'])>lastptid]

def saveTypeidfromAPI(savepath,nsearchmax=30,overwrite=True):

    if savepath==None:
        _savepath = getdatapath()+"typeid.ndjson"
    
    else:
        _savepath = savepath

    if overwrite:
        os.system("rm -rf "+_savepath)

    for nn in range(nsearchmax):
        typeid=nn+1
        url = "https://api.yamareco.com/api/v1/searchPoi"
        params = {"page":1,"type_id":typeid}
        res = requests.post(url,data=params)
        json_str = res.content.decode("euc-jp")
        json_dict = json.loads(json_str)

        try:
            typelist = json_dict['poilist'][0]['types']
            for i in range(len(typelist)):
                temp = typelist[i]
                temp.pop("detail")
                
                if (temp['type_id'] == str(typeid)):
                    
                    with open(_savepath,"a") as f:
                        writer = ndjson.writer(f)
                        writer.writerow(temp)
        except:
            break

def saveArealistfromAPI(savepath,overwrite=True):

    if savepath==None:
        _savepath = getdatapath()+"area.ndjson"
    
    else:
        _savepath = savepath

    if overwrite:
        os.system("rm -rf "+_savepath)

    url = "https://api.yamareco.com/api/v1/getArealist"
    res = requests.get(url)
    json_str = res.content.decode("euc-jp")
    json_dict = json.loads(json_str)

    with open(_savepath,"a") as f:
        writer = ndjson.writer(f)
        [writer.writerow(comp) for comp in json_dict['arealist']]

def saveGenrelistfromAPI(savepath,overwrite=True):

    if savepath==None:
        _savepath = getdatapath()+"genre.ndjson"
    
    else:
        _savepath = savepath

    if overwrite:
        os.system("rm -rf "+_savepath)

    url = "https://api.yamareco.com/api/v1/getGenrelist"
    res = requests.get(url)
    json_str = res.content.decode("euc-jp")
    json_dict = json.loads(json_str)

    with open(_savepath,"a") as f:
        writer = ndjson.writer(f)
        [writer.writerow(comp) for comp in json_dict['genrelist']]