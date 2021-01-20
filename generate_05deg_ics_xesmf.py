# import necessary packages
import xarray as xr
import numpy as np
import xesmf as xe
import os
import sys
import gridfill

# code is adpated from MOM6:MOM_state_initialization
def adjusthToFitBathymetry(ssh, h,bathy):
    hTolerance = 1.0e-10 #<  Tolerance to exceed adjustment criteria [Z ~> m]
    Angstrom_Z=1.0e-10
    nx=h.shape[2]
    ny=h.shape[1]
    nz=h.shape[0]
    eta=np.zeros([nz+1,ny,nx])
    eta[0,:,:]=ssh
    for k in range(nz):
        eta[k+1,:,:]=eta[k,:,:]-h[k,:,:]
    dz=eta[-1,:,:]+bathy
    for k in range(nz+1):
        eta[k,:,:]=np.where(-eta[k,:,:] > (bathy[:,:] + hTolerance),-bathy,eta[k,:,:])
    for k in range(nz-1,0,-1):
        h[k,:,:]=np.where(eta[k,:,:] < (eta[k+1,:,:] + Angstrom_Z),Angstrom_Z,eta[k,:,:] - eta[k+1,:,:])
    for i in range(nx):
        for j in range(ny):

    #   The whole column is dilated to accommodate deeper topography than
    # the bathymetry would indicate.
            if -eta[nz,j,i] < (bathy[j,i] - hTolerance):
                if eta[0,j,i] <= eta[nz,j,i]:
                    #for k in range(nz):
                    h[:,j,i] = (eta[1,j,i] + bathy[j,i]) / np.float(nz)
                else:
                    dilate = (eta[0,j,i] + bathy[j,i]) / (eta[0,j,i] - eta[nz,j,i])
                    #for k in range(nz):
                    h[:,j,i] = h[:,j,i] * dilate
    return h,eta

# main program   

# specify an output resolution
ires = 'mx025'
ores = "mx050"
# specify a date
nargs=len(sys.argv[:])
if (nargs != 2): 
   print('need to specify a date in YYYYMMDDHH format')
   os._exit(3)
cdate=sys.argv[1]
if (len(cdate)!=10):
   print('need to specify a date in YYYYMMDDHH format')
   os._exit(10)
print('processing ',cdate)

# specify a location to use
#nemsrc     = "/scratch2/BMC/gsienkf/Philip.Pegion/UFS-coupled/ICS/source/WeightGen/TTout/"
nemsrc     = "/work/noaa/marine/Jiande.Wang/UFS-python/FIX/WeightGen/TTout/"
# specifiy output directory
#outdir     = "/scratch2/BMC/gsienkf/Philip.Pegion/UFS-coupled/ICS/"+ores+"/"+cdate+"/"
outdir     = "/work/noaa/marine/Jiande.Wang/UFS-python/OUTPUT/"+ores+"/"+cdate+"/"

os.mkdir(outdir)
output_file= 'MOM6.mx050.ic.nc'
# ocean model restart location
#dirsrc     = "/scratch1/NCEPDEV/nems/Bin.Li/S2S/FROM_HPSS/"+cdate+"/mom6_da/"
#dirsrc     = "/work/noaa/marine/Jiande.Wang/UFS-python/SOURCE/"+cdate+"/mom6_da/"
dirsrc     = "/work/noaa/marine/Partha.Bhattacharjee/IC_Dir/CPC3Dvar/"+cdate+"/ocn/025/"

# target resolution bathymetry files
bathy_file='/work/noaa/marine/Jiande.Wang/UFS-python/FIX/topo/05deg_basedir/topog.nc'
#edits_file='/work/noaa/marine/Jiande.Wang/UFS-python/FIX/topo/1deg_basedir/topo_edits_011818.nc'

# OPEN ESMF weights file
wgt_dir='/work/noaa/marine/Jiande.Wang/UFS-python/FIX/WeightGen/TTout/'
grid_file_in=xr.open_dataset(wgt_dir+'tripole.mx025.nc')
grid_file_out=xr.open_dataset(wgt_dir+'tripole.mx050.nc')
wgtsfile_u_to_t = wgt_dir+'tripole.mx025.Cu.to.Ct.bilinear.nc'
wgtsfile_v_to_t = wgt_dir+'tripole.mx025.Cv.to.Ct.bilinear.nc'
wgtsfile_t_to_t = wgt_dir+'tripole.mx025.Ct.to.mx050.Ct.bilinear.nc'
wgtsfile_t_to_u = wgt_dir+'tripole.mx050.Ct.to.Cu.bilinear.nc'
wgtsfile_t_to_v = wgt_dir+'tripole.mx050.Ct.to.Cv.bilinear.nc'

# OPEN 1/4 degree initial condition files

res0=xr.open_dataset(dirsrc+"MOM.res.nc")
res1=xr.open_dataset(dirsrc+"MOM.res_1.nc")

# rename lat and lons from grid files for ESMF interpolation
mx025_t_grid=grid_file_in.rename({'lonCt': 'lon', 'latCt': 'lat'})
mx025_v_grid=grid_file_in.rename({'lonCv': 'lon', 'latCv': 'lat'})
mx025_u_grid=grid_file_in.rename({'lonCu': 'lon', 'latCu': 'lat'})

mx050_t_grid=grid_file_out.rename({'lonCt': 'lon', 'latCt': 'lat'})
mx050_v_grid=grid_file_out.rename({'lonCv': 'lon', 'latCt': 'lat'})
mx050_u_grid=grid_file_out.rename({'lonCu': 'lon', 'latCu': 'lat'})

# rotation angles for u and v currents
ang_in = grid_file_in.anglet
ang_out= grid_file_out.anglet

# open input data
res0=xr.open_dataset(dirsrc+"MOM.res.nc")
res1=xr.open_dataset(dirsrc+"MOM.res_1.nc")


rg_tt = xe.Regridder(mx025_t_grid, mx050_t_grid, 'bilinear',periodic=True,reuse_weights=True, filename=wgtsfile_t_to_t)
rg_ut = xe.Regridder(mx025_u_grid, mx025_t_grid, 'bilinear',periodic=True,reuse_weights=True, filename=wgtsfile_u_to_t)
rg_vt = xe.Regridder(mx025_v_grid, mx025_t_grid, 'bilinear',periodic=True,reuse_weights=True, filename=wgtsfile_v_to_t)
rg_tu = xe.Regridder(mx050_t_grid, mx050_u_grid, 'bilinear',periodic=True,reuse_weights=True, filename=wgtsfile_t_to_u)
rg_tv = xe.Regridder(mx050_t_grid, mx050_v_grid, 'bilinear',periodic=True,reuse_weights=True, filename=wgtsfile_t_to_v)

nx=len(grid_file_out.ni.values)
ny=len(grid_file_out.nj.values)
nz=len(res0.Layer.values)

# define land masks
lmask_in = xr.where(grid_file_in['wet'] > 0.0, 1.0, 0.0)
lmask_out = xr.where(grid_file_out['wet'] > 0.0, 1.0, 0.0)
# interpolate mask to new grid
lmask_interp0 = rg_tt(lmask_in.values)
lmask_interp = np.where(lmask_interp0 < 0.99, 0.0, 1.0)

#3d-copies
lmask_interp_3d=np.zeros([nz,ny,nx])
lmask_out_3d = np.zeros([nz,ny,nx])

for i in range(nz):
    lmask_interp_3d[i,:,:]=lmask_interp[:,:]
    lmask_out_3d[i,:,:]=grid_file_out['wet'].values

# interpolate values to new grid
new_t = rg_tt(res0['Temp'].values)
new_s = rg_tt(res0['Salt'].values)
new_h = rg_tt(res0['h'].values)
new_sfc = rg_tt(res1['sfc'].values)

# un-stagger currents
ut = rg_ut(res0['u'].values)
vt = rg_vt(res1['v'].values)

# rotate currents to earth relative
urot =   ut*np.cos(ang_in.values) +   vt*np.sin(ang_in.values)
vrot =   vt*np.cos(ang_in.values) -   ut*np.sin(ang_in.values)

# interpolate to new grid
new_urot = rg_tt(urot)
new_vrot = rg_tt(vrot)

# rotate currents back to grid relatvie (doesn't look right)
new_ut =   new_urot*np.cos(ang_out.values) - new_vrot*np.sin(ang_out.values)
new_vt =   new_vrot*np.cos(ang_out.values) + new_urot*np.sin(ang_out.values)

# re-stagger currents
new_u = rg_tu(new_ut)
new_v = rg_tv(new_vt)

# set options for poisson_grid_fill
guess     = True             # use zonal means
is_cyclic = True             # cyclic [global]
eps       = 1.e-2            # variable dependent
opt       = 0                # not used
relc      = 0.6              # relaxation coefficient
nscan     = 1500             # usually much less than this
# fill u v t etc.

ma=np.ma.masked_where(lmask_interp_3d == 0,new_t[0])
filled, converged = gridfill.fill(ma, 2,1,eps,relax=relc,itermax=nscan,initzonal=guess, cyclic=True)
filled = np.where(lmask_out_3d == 0.0, 0.0, filled)
new_t=filled.reshape([1,nz,ny,nx])

ma=np.ma.masked_where(lmask_interp_3d == 0,new_h[0])
filled, converged = gridfill.fill(ma, 2,1,eps,relax=relc,itermax=nscan,initzonal=guess, cyclic=True)
filled = np.where(lmask_out_3d == 0.0, 0.0, filled)
new_h=filled.reshape([1,nz,ny,nx])

ma=np.ma.masked_where(lmask_interp_3d == 0,new_s[0])
filled, converged = gridfill.fill(ma, 2,1,eps,relax=relc,itermax=nscan,initzonal=guess, cyclic=True)
filled = np.where(lmask_out_3d == 0.0, 0.0, filled)
new_s=filled.reshape([1,nz,ny,nx])

ma=np.ma.masked_where(lmask_interp_3d == 0,new_u[0])
filled, converged = gridfill.fill(ma, 2,1,eps,relax=relc,itermax=nscan,initzonal=guess, cyclic=True)
filled = np.where(lmask_out_3d == 0.0, 0.0, filled)
new_u=filled.reshape([1,nz,ny,nx])

ma=np.ma.masked_where(lmask_interp_3d == 0,new_v[0])
filled, converged = gridfill.fill(ma, 2,1,eps,relax=relc,itermax=nscan,initzonal=guess, cyclic=True)
filled = np.where(lmask_out_3d == 0.0, 0.0, filled)
new_v=filled.reshape([1,nz,ny,nx])

ma=np.ma.masked_where(lmask_interp == 0,new_sfc[0])
filled, converged = gridfill.fill(ma, 1,0,eps,relax=relc,itermax=nscan,initzonal=guess, cyclic=True)
filled = np.where(lmask_out == 0.0, 0.0, filled)
new_sfc=filled.reshape([1,ny,nx])

#open bathymerty file for vertical adjustment
bathy=xr.open_dataset(bathy_file)
#edits=xr.open_dataset(edits_file)
# edit bathymetry
#for i in edits.nEdits.values:
#    bathy.depth[edits.jEdit[i].values,edits.iEdit[i].values]=edits.zEdit[i].values
# set min/max of bathy
bathy.depth[:,:]=xr.where(bathy.depth < 9.5,9.5,bathy.depth)
bathy.depth[:,:]=xr.where(bathy.depth > 6500,6500,bathy.depth)
# adjust heights
new_h[0,:,:,:],eta=adjusthToFitBathymetry(new_sfc[0,:,:],new_h[0,:,:,:],bathy.depth.values)

# set up output file grid
new_lath=grid_file_out.latCt.mean(axis=1).values
new_lonh=grid_file_out.lonCt.mean(axis=0).values
new_latq=grid_file_out.latCv.mean(axis=1).values
new_lonq=grid_file_out.lonCu.mean(axis=0).values
# create xarray DataArrays
da_t = xr.DataArray(new_t,coords=({'lath' : (['lath'], new_lath), 'lonh' : (['lonh'], new_lonh), 'Layer': (['Layer'], res0['Layer'].values),\
                     'Time' : (['Time'], res0['Time'].values)}), dims=['Time','Layer','lath','lonh'])
da_s = xr.DataArray(new_s,coords=({'lath' : (['lath'], new_lath), 'lonh' : (['lonh'], new_lonh), 'Layer': (['Layer'], res0['Layer'].values), \
                     'Time' : (['Time'], res0['Time'].values)}), dims=['Time','Layer','lath','lonh'])
da_h = xr.DataArray(new_h,coords=({'lath' : (['lath'], new_lath), 'lonh' : (['lonh'], new_lonh), 'Layer': (['Layer'], res0['Layer'].values), \
                     'Time' : (['Time'], res0['Time'].values)}), dims=['Time','Layer','lath','lonh'])
da_sfc=xr.DataArray(new_sfc,coords=({'lath' : (['lath'], new_lath), 'lonh' : (['lonh'], new_lonh), 'Time' : (['Time'], res0['Time'].values)}), \
                     dims=['Time','lath','lonh'])
da_u = xr.DataArray(new_u,coords=({'lath' : (['lath'], new_lath), 'lonq' : (['lonq'], new_lonq), 'Layer': (['Layer'], res0['Layer'].values), \
                     'Time' : (['Time'], res0['Time'].values)}), dims=['Time','Layer','lath','lonq'])
da_v = xr.DataArray(new_v,coords=({'latq' : (['latq'], new_latq), 'lonh' : (['lonh'], new_lonh), 'Layer': (['Layer'], res0['Layer'].values), \
                     'Time' : (['Time'], res0['Time'].values)}), dims=['Time','Layer','latq','lonh'])
# create xarray DataSets
ds_t=da_t.to_dataset(name='Temp')
ds_s=da_s.to_dataset(name='Salt')
ds_h=da_h.to_dataset(name='h')
ds_sfc=da_sfc.to_dataset(name='sfc')
ds_u=da_u.to_dataset(name='u')
ds_v=da_v.to_dataset(name='v')
# merge variables
ds_out=xr.merge([ds_t,ds_s,ds_h,ds_sfc,ds_u,ds_v])
# write to file
ds_out.to_netcdf(outdir+output_file,unlimited_dims='Time')
