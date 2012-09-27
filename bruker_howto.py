# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 13:03:52 2012

Import Bruker howto script.
Sorry, there really isn't much here, but this is to just give the general outline
of how to go about reconstructing Bruker data.

@author: eric
"""

from read_bruker import read_bruker, read_bruker_text, read_bruker_binary, grid_bruker
from os.path import join, dirname, realpath
from matplotlib.pyplot import show, imshow, figure

scan_directory='some path'
params={}
params['scan']=read_bruker_text(join(scan_directory,'method')) #get information about the scan
matrix=[int(params['scan']['PVM_Matrix'][0]),int(params['scan']['PropellerBladeRes']),len(params['scan']['PVM_EncChanScaling']),int(params['scan']['PVM_SPackArrNSlices'][0]),int(params['scan']['PropellerNbAngleInc'])]
params['fid']=read_bruker_binary(join(scan_directory,'fid'),matrix,'complex',1) #get the raw data
img_square=grid_bruker(params,'bladewise') #IT IS UP TO YOU TO RECONSTRUCT AS YOU WANT!

figure()
imshow(abs(img_square))
show()
