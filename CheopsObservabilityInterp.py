from astropy.coordinates import SkyCoord,GeocentricMeanEcliptic
from astropy import units as u
from scipy.interpolate import CloughTocher2DInterpolator as ct2d
import numpy as np
import pandas as pd


def InterpCheopsObs(radec,obs_interps=None):
    # Uses pre-computed interpolated function of the Cheops Sky (https://www.cosmos.esa.int/web/cheops/the-cheops-sky)
    # to estimate Cheops observability given position
    # 
    # INPUTS:
    #  SkyCoord ra/dec stellar coordinates
    #  List of interpolation functions (output of init_cheops_func() - necessary if looping over many stars)
    # 
    # OUTPUTS:
    #  cumulative time observable at 50% (e.g. N_days_observable*efficiency)
    #  total number of days the object is observable at 50%
    #  Max efficiency (50, 70 or 80%)
    
    if obs_interps is None:
        obs_interps=init_cheops_func()
    
    sc=radec.transform_to(GeocentricMeanEcliptic)
    lon=sc.lon.deg;lat=sc.lat.deg
    obsty=np.array([obs_interps[0](lon,lat),obs_interps[1](lon,lat),obs_interps[2](lon,lat)])

    #We need to convert from "cumulative time" to "N_days observable"
    # We do this by summing observable days at >80% (ssuming 90% eff), at 70-80 (assume 75%) & at 50% (assume 60%)
    N_days_observable = np.clip((obsty[0]-obsty[1])/0.60 + (obsty[1]-obsty[2])/0.75 + obsty[2]/0.90, 0.0, 365)
    #print(obsty,n_days_obs)
    
    max_eff=80 if obsty[2]>0.5 else 70 if obsty[1]>0.5 else 50 if obsty[0]>0.5 else 0.0
    
    return obsty[0], N_days_observable, max_eff

def init_cheops_func():

    #50% efficiency data:
    cheopsRA = np.loadtxt("data/ra_grid.dat")
    cheopsDEC = np.loadtxt("data/dec_grid.dat")
    cheopsSky = np.loadtxt("data/6am_700_10_noMoon_conf5_V9_50d_49m.dat")

    cheopsSky_nonan=np.nan_to_num(cheopsSky,0.0)

    scs_50=SkyCoord((cheopsRA.ravel()%360)*u.deg,cheopsDEC.ravel()*u.deg)
    scs2_50=scs_50.transform_to(GeocentricMeanEcliptic)
    xyz_50=np.column_stack((scs2_50.lon.deg,scs2_50.lat.deg,cheopsSky_nonan.ravel()))
    mirror_neg=xyz_50[:,0]>350
    mirror_pos=xyz_50[:,0]<10
    xyz_50=np.vstack((np.column_stack((xyz_50[mirror_neg,0]-360,xyz_50[mirror_neg,1:])),
                      xyz_50,
                      np.column_stack((xyz_50[mirror_pos,0]+360,xyz_50[mirror_pos,1:]))))
    obs_interp_50i=ct2d(xyz_50[:,:2],xyz_50[:,2])
    obs_interp_50=lambda lon,lat: np.nan_to_num(np.clip(obs_interp_50i(lon,lat),0.0,80),0.0)
    
    #70% efficiency data:
    data_70 = pd.read_csv("data/70pc_Cheops_Map_data.csv",header=None).values
    data_70[:,0]=(-15*data_70[:,0])%360

    scs_70=SkyCoord(data_70[:,0]*u.deg,data_70[:,1]*u.deg)
    scs2_70=scs_70.transform_to(GeocentricMeanEcliptic)

    xyz_70=np.column_stack((scs2_70.lon.deg,scs2_70.lat.deg,data_70[:,2]))
    mirror_neg=xyz_70[:,0]>350
    mirror_pos=xyz_70[:,0]<10
    xyz_70=np.vstack((np.column_stack((xyz_70[mirror_neg,0]-360,xyz_70[mirror_neg,1:])),
                      xyz_70,
                      np.column_stack((xyz_70[mirror_pos,0]+360,xyz_70[mirror_pos,1:]))))
    obs_interp_70i=ct2d(xyz_70[:,:2],xyz_70[:,2])
    obs_interp_70=lambda lon,lat: np.nan_to_num(np.clip(obs_interp_70i(lon,lat),0.0,80),0.0)
    
    #80% efficiency data:
    data_80 = pd.read_csv("data/80pc_Cheops_Map_data.csv",header=None).values
    data_80[:,0]=(-15*data_80[:,0])%360

    scs_80=SkyCoord(data_80[:,0]*u.deg,data_80[:,1]*u.deg)
    scs2_80=scs_80.transform_to(GeocentricMeanEcliptic)

    xyz_80=np.column_stack((scs2_80.lon.deg,scs2_80.lat.deg,data_80[:,2]))
    mirror_neg=xyz_80[:,0]>350
    mirror_pos=xyz_80[:,0]<10
    xyz_80=np.vstack((np.column_stack((xyz_80[mirror_neg,0]-360,xyz_80[mirror_neg,1:])),
                      xyz_80,
                      np.column_stack((xyz_80[mirror_pos,0]+360,xyz_80[mirror_pos,1:]))))

    obs_interp_80i=ct2d(xyz_80[:,:2],xyz_80[:,2])
    obs_interp_80=lambda lon,lat: np.nan_to_num(np.clip(obs_interp_80i(lon,lat),0.0,80),0.0)
    return [obs_interp_50,obs_interp_70,obs_interp_80]
