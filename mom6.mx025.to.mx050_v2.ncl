;------------------------------------------------------------------
; Denise.Worthen@noaa.gov (Feb 2019)

  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;----------------------------------------------------------------------
begin

   ; specify a date
            cdate = "2016100100"
   ; specify a location to use
       nemsrc     = "/scratch2/NCEPDEV/climate/Denise.Worthen/TTout2/"
   ; specifiy output directory 
       outdir     = "/scratch2/NCEPDEV/climate/Denise.Worthen/ICmx050/"+cdate+"/"
   ; interpolation methods
   methods        = (/"bilinear"/)
   ; ocean model restart location 
       dirsrc     = "/scratch1/NCEPDEV/nems/Bin.Li/S2S/FROM_HPSS/"+cdate+"/mom6_da/"

   restartfilelist = systemfunc("ls "+dirsrc+"MOM.res*.nc")
          restartf = addfiles(restartfilelist,"r")
            nfiles = dimsizes(restartfilelist)
   ;print(restartfilelist)
 
   gridsrc = nemsrc+"tripole.mx025.nc"
   ; create a temporary restart file with needed variables at
   ; qdeg resolution
   outfile = outdir+"MOM6.mx025.res.nc"

   if(.not.fileexists(outfile))then
     print("creating intermediate restart file "+outfile)
     ncocmd = "ncks -O -v Temp,Salt,h,u "+restartfilelist(0)
     print(ncocmd+"  "+outfile)
    system(ncocmd+"  "+outfile)
    ; v-velocity in _1 file
     ncocmd = "ncks -A -v v,sfc "+restartfilelist(1)
     print(ncocmd+"  "+outfile)
    system(ncocmd+"  "+outfile)
   else
     print("intermediate restart file exists, continuing ")
   end if
   delete(restartf)
   delete(gridsrc)

    ; set options for poisson_grid_fill
    guess     = 1                ; use zonal means
    is_cyclic = True             ; cyclic [global]
    nscan     = 1500             ; usually much less than this
    eps       = 1.e-2            ; variable dependent
    relc      = 0.6              ; relaxation coefficient
    opt       = 0                ; not used

;----------------------------------------------------------------------
;
;----------------------------------------------------------------------

   restartf = addfile(outfile,"r")
   ;print(restartf)
   
   vnames = getfilevarnames(restartf)
    nvars = dimsizes(vnames)
   dNames = getfiledimnames(restartf)
   dSizes = getfiledimsizes(restartf)
   ;print("file "+outfile+" contains "+dNames+"    "+dSizes)

   ; load the time and layer values, these don't change
     Time = restartf->Time
    Layer = restartf->Layer
    
   ; output resolution sizes is known
      nt = 1
      nz = dSizes(0)
      nj = 576
      ni = 720

   gridtype = new((/nvars/),"string")
   gridtype = ""
    nvpairs = 1
    veclist = new( (/nvpairs,2,2/),"string")
    veclist = (/ (/  "u",   "v"/), (/"Cu", "Cv"/) /) 
    print(veclist)

   ; use the land mask to set the masking
   ; the mask contain 1's on land and 0's at valid points.
       gf = addfile(nemsrc+"tripole.mx025.nc","r")
    tmask = where(gf->wet .eq. 0.0, 1.0, 0.0)
   mask2d = conform_dims(dimsizes(restartf->sfc), tmask, (/1,2/))
   mask3d = conform_dims(dimsizes(restartf->h),   tmask, (/2,3/))

   ; use the interface thickness to set the bathy masking
   ; the mask contain 1's on land and 0's at valid points.
     minh = 1.0e-3       ; the minimum thickness allowed in MOM6
    bathy = where(restartf->h .le. minh, 1.0, 0.0)
   delete(tmask)
   delete(gf)

   ;printVarSummary(mask2d)
   ;printVarSummary(mask3d)
   ; regrid the masks to obtain the interpolation masks.
   ; the mask2d,mask3d contain 1's on land and 0's at valid points.
   ; when remapped, any mask value > 0 identifies land values that
   ; have crept into the field. remapped model fields are then
   ; masked with this interpolation mask
   ; specify the tripole-to-tripole weights file
   wgtsfile = nemsrc+"tripole.mx025.Ct.to.mx050.Ct.bilinear.nc"
   rgmask2d = ESMF_regrid_with_weights(mask2d, wgtsfile,False)
   rgmask3d = ESMF_regrid_with_weights(mask3d, wgtsfile,False)
   rgbath3d = ESMF_regrid_with_weights(bathy, wgtsfile,False)

    ; create the interpolation mask
    rgmask2d = where(rgmask2d .gt. 0.0, default_fillvalue(typeof(rgmask2d)), 1.0)
    rgmask3d = where(rgmask3d .gt. 0.0, default_fillvalue(typeof(rgmask3d)), 1.0)
    rgbath3d = where(rgbath3d .gt. 0.0, default_fillvalue(typeof(rgbath3d)), 1.0)

   ; retrieve the rotation angle on the source grid and conform to 3d
         gf = addfile(nemsrc+"tripole.mx025.nc","r")
     angsrc = gf->anglet
   angsrc3d = conform_dims(dimsizes(mask3d),angsrc,(/2,3/))
   delete(gf)
   ;printVarSummary(angsrc3d)

   ; retrieve the rotation angle on the destination grid and conform to 3d
         gf = addfile(nemsrc+"tripole.mx050.nc","r")
     angdst = gf->anglet
   angdst3d = conform_dims(dimsizes(rgmask3d),angdst,(/2,3/))

   ; retrieve known land mask on destination grid and conform to 3d
    tmask = gf->wet
   land2d = conform_dims((/nt,nj,ni/),    tmask, (/1,2/)) 
   land3d = conform_dims((/nt,nz,nj,ni/), tmask, (/2,3/))
   ;printVarSummary(land3d)
   delete(tmask)
   delete(gf)

    testfile = "mask_bilinear.nc"
    system("/bin/rm -f "+testfile)
    ; create
    testcdf = addfile(testfile,"c")
    testcdf->rgmask2d = rgmask2d
    testcdf->rgmask3d = rgmask3d
    testcdf->mask3d = mask3d
    testcdf->bathy = bathy
    testcdf->land3d = land3d
    testcdf->rgbath3d = rgbath3d
    ; close
    ;delete(testcdf)

;----------------------------------------------------------------------
; a set of values for the dimension variables on the 050 grid
;----------------------------------------------------------------------

      gf = addfile(nemsrc+"tripole.mx050.nc","r")
   latCt = gf->latCt(:,ni/4)
   lonCt = where(gf->lonCt(ni/2,:) .lt. 60.0d0, gf->lonCt(ni/2,:), gf->lonCt(ni/2,:)-360.d0)
   latBu = gf->latBu(:,ni/4)
   lonBu = where(gf->lonBu(ni/2,:) .lt. 60.5d0, gf->lonBu(ni/2,:), gf->lonBu(ni/2,:)-360.d0)
   delete(gf)

   copy_VarAtts(restartf->lath,latCt) 
   copy_VarAtts(restartf->lonh,lonCt) 
   copy_VarAtts(restartf->latq,latBu) 
   copy_VarAtts(restartf->lonq,lonBu) 
   ;printVarSummary(latCt)
   ;printVarSummary(lonCt)
   ;printVarSummary(latBu)
   ;printVarSummary(lonBu)

;   latCt = ispan(1,nj,1)
;   lonCt = ispan(1,ni,1)
;   latBu = ispan(1,nj,1)
;   lonBu = ispan(1,ni,1)
;   printVarSummary(latCt)

;----------------------------------------------------------------------
; set up the output netcdf file
;----------------------------------------------------------------------

    ;rgrestartf = nemsrc+"MOM6.mx050.res.nc"
    rgrestartf = outdir+"MOM6.mx050.res.nc"
    print("regridded output file name "+rgrestartf)
    system("/bin/rm -f " + rgrestartf)    ; remove if exists
    outcdf  = addfile (rgrestartf, "c")  ; open output file
    ;print(outcdf)

    ; explicitly declare file definition mode. Improve efficiency.
    setfileoption(outcdf,"DefineMode",True)

    ; create global attributes of the file
    fAtt               = True            ; assign file attributes
    fAtt@creation_date = systemfunc ("date")
    fAtt@source_file   = outfile
    fileattdef( outcdf, fAtt )           ; copy file attributes

    ; predefine the coordinate variables and their dimensionality
    dimNames = dNames
    ;print(dimNames)
    dimNamesReorder = (/dNames(1),dNames(0),dNames(2),dNames(3),dNames(4),dNames(5)/)
    ;print(dimNamesReorder)

    dimSizes = (/   1,     nz,    nj,     ni,    ni,    nj/)
    dimUnlim = (/True,  False, False,  False, False, False/)
    filedimdef(outcdf,dimNamesReorder,dimSizes,dimUnlim)
    ;print(dimSizes)

    ; predefine each variable in the output file from the original restart
    do nv = 0,nvars-1
      varname = vnames(nv)
       vtype := getfilevartypes(restartf,varname)
      vatts  := getfilevaratts(restartf,varname)       ; The := is important because these
      dnames := getfilevardimnames(restartf,varname)   ; variables may be a different size
      dsizes := getfilevardimsizes(restartf,varname)   ; the next time through the loop
       ndims := dimsizes(dnames)
      ;print(nv+"  "+varname+" "+ndims+" "+dnames)

      ; set grid stagger for 3d fields
      if(ndims .ge. 3)then
       ;print(varname+" "+ndims+" "+dnames(ndims-2)+" "+dnames(ndims-1))
       if(dnames(ndims-2) .eq. "lath" .and. dnames(ndims-1) .eq. "lonh")then
        gridtype(nv) = "Ct"
       end if
       if(dnames(ndims-2) .eq. "lath" .and. dnames(ndims-1) .eq. "lonq")then
        gridtype(nv) = "Cu"
       end if
       if(dnames(ndims-2) .eq. "latq" .and. dnames(ndims-1) .eq. "lonh")then
        gridtype(nv) = "Cv"
       end if
       ;print(varname+" "+ndims+"  "+gridtype(nv))
      end if

      filevardef(outcdf, varname, vtype, dnames)
      ;print(varname+"  "+gridtype(nv))
    end do
    ; explicitly exit file definition mode.
    setfileoption(outcdf,"DefineMode",False)

    ; write the dimensions to the file
    outcdf->Time    = (/Time/)
    outcdf->Layer   = (/Layer/)
    ; and the dimension arrays
    outcdf->lonh    = (/lonCt/)
    outcdf->lath    = (/latCt/)
    outcdf->lonq    = (/lonBu/)
    outcdf->latq    = (/latBu/)

;----------------------------------------------------------------------
; loop over nvars variables, doing only 3d (sfc) variables
;----------------------------------------------------------------------

    do nv = 0,nvars-1
     varname = vnames(nv)
     vargrid = gridtype(nv)
      ocnvar = restartf->$varname$
       ndims = dimsizes(dimsizes(ocnvar))
     if(ndims .eq. 3 .and. vargrid .eq. "Ct")then
      print("regridding variable "+varname+"  "+ndims)

      rgtt = ESMF_regrid_with_weights(ocnvar,wgtsfile,False)
      delete(rgtt@lon2d)
      delete(rgtt@lat2d)
      ; enter file definition mode to add variable attributes
      setfileoption(outcdf,"DefineMode",True)
      filevarattdef(outcdf, varname, rgtt)
      setfileoption(outcdf,"DefineMode",False)

      ; mask out land contaminated areas
      rgtt = where(ismissing(rgmask2d), rgtt@_FillValue, rgtt)
      ; fill missing areas
      print("filling variable "+varname)
      poisson_grid_fill( rgtt, is_cyclic, guess, nscan, eps, relc, opt)

      ; mask with land mask
      rgtt = where(land2d .eq. 0.0, 0.0, rgtt)
      outcdf->$varname$  = (/rgtt/)

      delete(rgtt)
     end if
     delete(ocnvar)
    end do

;----------------------------------------------------------------------
; loop over nvars variables, not vectors or 1-dims
;----------------------------------------------------------------------

    ;nv = 4; h
    do nv = 0,nvars-1
     varname = vnames(nv)
     vargrid = gridtype(nv)
      ocnvar = restartf->$varname$
       ndims = dimsizes(dimsizes(ocnvar))
     if(ndims .ge. 4 .and. vargrid .eq. "Ct")then
      print("regridding variable "+varname+"  "+ndims)

      rgtt = ESMF_regrid_with_weights(ocnvar,wgtsfile,False)
      delete(rgtt@lon2d)
      delete(rgtt@lat2d)
      ; enter file definition mode to add variable attributes
      setfileoption(outcdf,"DefineMode",True)
      filevarattdef(outcdf, varname, rgtt)
      setfileoption(outcdf,"DefineMode",False)

      ; mask out land contaminated areas
       rgtt = where(ismissing(rgmask3d), rgtt@_FillValue, rgtt)
      ; fill missing areas
      poisson_grid_fill(rgtt, is_cyclic, guess, nscan, eps, relc, opt)

      ; mask with land mask
      ; if this is thickness, set min val
      if(varname .eq. "h")then
       rgtt = where(ismissing(rgbath3d) .and. land3d .eq. 0.0, minh, rgtt)
      else
       rgtt = where(land3d .eq. 0.0, 0.0, rgtt)
      end if 
      outcdf->$varname$  = (/rgtt/)

      delete(rgtt)
     end if
     delete(ocnvar)
    end do

;----------------------------------------------------------------------
; retrieve and regrid the vectors
;----------------------------------------------------------------------

   do nv = 0,nvpairs-1
     vecnames = veclist(nv,0,:)
     vecgrids = veclist(nv,1,:)

     ; unstagger velocities from (Cu,Cv) to Ct
            n = 0
         uvel = restartf->$vecnames(n)$
       vecfld = where(ismissing(uvel),0.0,uvel)
       copy_VarAtts(uvel,vecfld)
     print("unstagger "+vecnames(n)+" from "+vecgrids(n)+" to Ct")
     wgtsfile = nemsrc+"tripole.mx025."+vecgrids(n)+".to.Ct.bilinear.nc"
           ut = ESMF_regrid_with_weights(vecfld,wgtsfile,False)

            n = 1
         vvel = restartf->$vecnames(n)$
       vecfld = where(ismissing(vvel),0.0,vvel)
       copy_VarAtts(vvel,vecfld)
     print("unstagger "+vecnames(n)+" from "+vecgrids(n)+" to Ct")
     wgtsfile = nemsrc+"tripole.mx025."+vecgrids(n)+".to.Ct.bilinear.nc"
           vt = ESMF_regrid_with_weights(vecfld,wgtsfile,False)
     delete(vecfld)

     ; rotate IJ to EW
     ; first copy Metadata
       urot = ut
       vrot = vt
       urot =   ut*cos(angsrc3d) +   vt*sin(angsrc3d)
       vrot =   vt*cos(angsrc3d) -   ut*sin(angsrc3d)
     ; copy back
        ut = urot
        vt = vrot
     delete([/urot, vrot/])
     ;printVarSummary(ut)

     ; regrid EW velocities on Ct grid to mx050 
     wgtsfile = nemsrc+"tripole.mx025.Ct.to.mx050.Ct.bilinear.nc"
            n = 0
       vecfld = where(ismissing(ut),0.0,ut)
      print("remapping EW velocity "+vecnames(n)+" on Ct")
      rgu = ESMF_regrid_with_weights(vecfld,wgtsfile,False)
      ; mask out land contaminated areas
      rgu = where(ismissing(rgmask3d), rgu@_FillValue, rgu)
      ; fill missing areas
      print("filling variable "+vecnames(n))
      poisson_grid_fill( rgu, is_cyclic, guess, nscan, eps, relc, opt)
      ; mask with land mask
      rgu = where(land3d .eq. 0.0, 0.0, rgu)

            n = 1
       vecfld = where(ismissing(vt),0.0,vt)
      print("remapping EW velocity "+vecnames(n)+" on Ct")
      rgv = ESMF_regrid_with_weights(vecfld,wgtsfile,False)
      ; mask out land contaminated areas
      rgv = where(ismissing(rgmask3d), rgv@_FillValue, rgv)
      ; fill missing areas
      print("filling variable "+vecnames(n))
      poisson_grid_fill( rgv, is_cyclic, guess, nscan, eps, relc, opt)
      ; mask with land mask
      rgv = where(land3d .eq. 0.0, 0.0, rgv)
      delete(vecfld)

     ; rotate EW to IJ
     ; first copy Metadata
       urot = rgu
       vrot = rgv
       urot =   rgu*cos(angdst3d) -   rgv*sin(angdst3d)
       vrot =   rgv*cos(angdst3d) +   rgu*sin(angdst3d)
     ; copy back
        rgu = urot
        rgv = vrot
     delete([/urot, vrot/])

     ; restagger from Ct to (Cu,Cv) 
            n = 0
       vecfld = where(ismissing(rgu),0.0,rgu)
       copy_VarAtts(rgu,vecfld)
     print("restagger "+vecnames(n)+" from Ct to "+vecgrids(n))
     wgtsfile = nemsrc+"tripole.mx050.Ct.to."+vecgrids(n)+".bilinear.nc"
            u = ESMF_regrid_with_weights(vecfld,wgtsfile,False)
     delete(u@lon2d)
     delete(u@lat2d)
     ; do not mask with land mask because land is on Ct
     ;u = where(land3d .eq. 0.0, 0.0, u)

     ; enter file definition mode to add variable attributes
     setfileoption(outcdf,"DefineMode",True)
     filevarattdef(outcdf, vecnames(n),   u)
     setfileoption(outcdf,"DefineMode",False)
     outcdf->$vecnames(n)$   = (/u/)

            n = 1
       vecfld = where(ismissing(rgv),0.0,rgv)
       copy_VarAtts(rgv,vecfld)
     print("restagger "+vecnames(n)+" from Ct to "+vecgrids(n))
     wgtsfile = nemsrc+"tripole.mx050.Ct.to."+vecgrids(n)+".bilinear.nc"
            v = ESMF_regrid_with_weights(vecfld,wgtsfile,False)
     delete(v@lon2d)
     delete(v@lat2d)
     ; do not mask with land mask because land is on Ct
     ;v = where(land3d .eq. 0.0, 0.0, v)

     ; enter file definition mode to add variable attributes
     setfileoption(outcdf,"DefineMode",True)
     filevarattdef(outcdf, vecnames(n),   v)
     setfileoption(outcdf,"DefineMode",False)
     outcdf->$vecnames(n)$   = (/v/)

     end do
exit
end
