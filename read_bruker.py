# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 13:55:50 2012

read_bruker

This file reads and reconstructs bruker propeller data.

@author: eric
"""

from matplotlib.pyplot import imshow, gca, show, figure, plot, subplot, axis, figtext
from numpy import zeros, exp, swapaxes, shape, fromfile, int32, reshape, prod, append, linspace, pi, array
from numpy.fft import fft, fftshift
from copy import copy
from phantom import phantom
from csv import reader
from os.path import join, dirname, realpath
import sys
from nfft_wrappers import init_nfft_2d, init_iterative_nfft_2d, grid_nfft_2d, iterative_nfft_2d, finalize_nfft_2d, finalize_iterative_nfft_2d
from nfft_helpers import propeller_sampling_location_generator, calc_density_weights, comp_fov


def read_bruker(directory):
    params={} #create dictionary for bruker files
    #b['reco']={}
    params['reco']=read_bruker_text(join(directory,'reco'))
    params['method']=read_bruker_text(join(dirname(dirname(directory)),'method'))
    fid=open(join(directory,'2dseq'),'rb')
    params['img']=fromfile(file=fid,dtype=int32)
    params['img']=reshape(params['img'],[int(params['reco']['RECO_size'][0]),int(params['reco']['RECO_size'][1]),int(params['reco']['RECO_size'][2]),len(params['method']['PVM_DwBvalEach'])+1],'F')
    return params
    

    
    
def read_bruker_text(filename):
    r=reader(open(filename,'r'),delimiter=' ')
    #b['reco']={}
    b={}
    for row in r:
        #print row
        if row[0][:3]=='##$': #name of parameter
            rs=row[0][3:].split('=')
            #print rs
            if rs[1]=='(': #if we have an array following
                try: #try to convert to floats, if possible
                    b[last_parameter]=[float(x) for x in b[last_parameter]]
                except:
                    pass
                last_parameter=rs[0]
            else:
                try:
                    b[rs[0]]=float(rs[1])
                except:
                    b[rs[0]]=rs[1]
        elif row[0][:2]=='##': #user info
            rs=row[0][2:].split('=')
            b[rs[0]]=rs[1]
        elif row[0][:2]!='$$': #as long as it isn't a comment
            #print row
            #print last_parameter
            try:
                b[last_parameter].extend(filter(None,row)) #filter to prevent empties
            except:
                b[last_parameter]=filter(None,row)
            #print row
            #print row[0][3:].split('=')
        #print row[0][:3]
    return b
        
        
def read_bruker_2dseq(directory):
    b={} #create dictionary for bruker files
    #b['reco']={}
    b['reco']=read_bruker_text(join(directory,'reco'))
    b['method']=read_bruker_text(join(dirname(dirname(directory)),'method'))
    fid=open(join(directory,'2dseq'),'rb')
    b['img']=fromfile(file=fid,dtype=int32)
    b['img']=reshape(b['img'],[int(b['reco']['RECO_size'][0]),int(b['reco']['RECO_size'][1]),int(b['reco']['RECO_size'][2]),len(b['method']['PVM_DwBvalEach'])+1],'F')
    return b
    
    
def read_bruker_binary(filename,matrix,complex_flag='real',flipodd=0):
    #returns: x,y,blade,coil,slice
    fid=open(filename,'rb')
    if complex_flag=='real':
        img=reshape(fromfile(file=fid,dtype=int32),matrix,'F')
    elif complex_flag=='complex':
        #img=zeros(matrix,complex)
        timg=fromfile(file=fid,dtype=int32)
        #img=timg[0::2].astype(float)+1j*timg[1::2].astype(float)
        img=timg[0::2]+1j*timg[1::2]
        if (matrix[0]*matrix[1]*2 % 1024): #some strange bruker 1024 byte thing
            print('CAREFUL! The blade is not a multiple of 1024, data may not be correct!')
        print(shape(timg), prod(matrix)*2)
        img=reshape(img,matrix,'F')
    if flipodd==1:
        img[:,::2,:,:,:]=copy(img[::-1,::2,:,:,:])
    return img
    

    
    
    
def grid_bruker(dic,option='bladewise'):
    nangles=int(dic['scan']['PropellerNbAngleInc'])
    blade=map(int,[dic['scan']['PVM_EncMatrix'][0],dic['scan']['PropellerBladeRes']])
    #print type(blade)
    #print blade
    #print nangles
    
    fov=comp_fov(dic['scan']['PVM_EncMatrix'][0])
    #print shape(locs)
    fid=dic['fid'].squeeze() #just squeeze out other dims for now!
    
    try: #try to use external calibration data
        angles=linspace(0,2*pi/dic['scan']['PropellerUndersampling'],nangles,False)
        #print angles
        traj_x=array(dic['scan']['external_trajectory']['PVM_EpiTrajAdjkx'])/blade[0]
        traj_y=array(dic['scan']['external_trajectory']['PVM_EpiTrajAdjky'])/blade[0]
        traj_x-=traj_x[blade[0]/2]
        traj_y-=traj_y[blade[0]/2]
        locs_vec,locs=propeller_sampling_location_generator(blade,nangles,x_true=traj_x,y_true=traj_y)
        #print traj_x
    except:
        locs_vec,locs=propeller_sampling_location_generator(blade,nangles)
    locs=reshape(locs,shape(fid),'F')
    #print traj_y
    if option=='bladewise':
        img_square=zeros(append(dic['scan']['PVM_EncMatrix'],nangles),complex)
        for ang in range(nangles):
            #idxs=range(ang*blade[1],(ang+1)*blade[1],1)
            #print idxs
            p=init_nfft_2d(locs[:,:,ang],int(dic['scan']['PVM_EncMatrix'][0]),fov)
            img_square[:,:,ang]=grid_nfft_2d(p,fid[:,:,ang])
            finalize_nfft_2d(p)
    elif option=='bladecombine':
        p=init_nfft_2d(locs,int(dic['scan']['PVM_EncMatrix'][0]),fov)
        img_square=grid_nfft_2d(p,dic['fid'])
        finalize_nfft_2d(p)
    return img_square
    