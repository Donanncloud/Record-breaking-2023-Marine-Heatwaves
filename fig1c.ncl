    latS   =  -60
    latN   =60
    lonL   =  0
    lonR   =  360

    ; northatlantic
   latS1   =  5
   latN1   =45
   lonL1   =  300
   lonR1   =  350

   ;tropical pacific

   latS2   =  -15
   latN2   =  15
   lonL2   =  210
   lonR2   =  280



   ;North pacific
   latS3   =  20
   latN3   =  60
   lonL3   =  120
   lonR3   =  240


   ;Southwest Pacific
   latS4   =  -58
   latN4   =  -38
   lonL4   =  150
   lonR4   =  220



    ;================================================
    diri  = "/public/home/dongty/code/2023global/"
    fili1 = "OISST_4mhw_counts_durs_mean-cum_global_allyear_90th_1982to2023.nc"  
    ; fili1="OISST_1x1_4mhw_counts_durs_mean-cum_global_allyear_90th_1982to2023.nc" 


 

    a1    = addfile(diri+fili1,"r")


    ccl     = a1->ccl(3,:,{latS:latN},{lonL:lonR})


    lat=a1->lat({latS:latN})
    lon=a1->lon({lonL:lonR})

    lat1=a1->lat({latS1:latN1})
   lon1=a1->lon({lonL1:lonR1})


   lat2=a1->lat({latS2:latN2})
   lon2=a1->lon({lonL2:lonR2})


   lat3=a1->lat({latS3:latN3})
   lon3=a1->lon({lonL3:lonR3})


   lat4=a1->lat({latS4:latN4})
   lon4=a1->lon({lonL4:lonR4})


    nlat=dimsizes(lat)
    nlon=dimsizes(lon)
    lat@units="degrees-north"
    lon@units="degrees-east"

    year=ispan(1982,2023,1)
    nyear=42



    ccl!0="year"
    ccl!1="lat"
    ccl!2="lon"
    ccl&year=year
    ccl&lat=lat
    ccl&lon=lon
    ccl&lat@units="degrees-north"
    ccl&lon@units="degrees-east"
    printVarSummary(ccl)


    ;space 
    ccl_sp=ccl(41,:,:)

    ccl_sp:=where(ccl_sp.ge.40.0,ccl_sp,ccl_sp@_FillValue)

    ccl_sp!0="lat"
    ccl_sp!1="lon"
    ccl_sp&lat=lat
    ccl_sp&lon=lon
    ccl_sp&lat@units="degrees-north"
    ccl_sp&lon@units="degrees-east"


;累加
   re   = 6.37122e03
   rad  = 4.0 * atan(1.0) / 180.0
   con  = re * rad                 
   clat = cos(lat * rad)           ; cosine of latitude
   
   dlon = (lon(2) - lon(1))        ; assume dlon is constant
   dlat = (lat(2) - lat(1))        ; assume dlat is constant

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dxdy = dy * dx                  ; dydx(nlat)
   

   
   wgt  = new((/nlat, nlon/), typeof(ccl_sp))
   wgt  = conform(wgt, dxdy, 0)
    
   ccl_global =wgt_areasum2(ccl_sp, wgt, 0)


    print(ccl_global)

    delete([/clat,dx,dy,dxdy,wgt/])

    ;1
   clat = cos(lat1 * rad)           ; cosine of latitude

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dxdy = dy * dx                  ; dydx(nlat)
   
   nlat1=dimsizes(lat1)
   nlon1=dimsizes(lon1)
   
   wgt  = new((/nlat1, nlon1/), typeof(ccl_sp({latS1:latN1},{lonL1:lonR1})))
   wgt  = conform(wgt, dxdy, 0)
    
   ccl_1 =wgt_areasum2(ccl_sp({latS1:latN1},{lonL1:lonR1}), wgt, 0)

   print(ccl_1)
   delete([/clat,dx,dy,dxdy,wgt/])

    ;2
   clat = cos(lat2 * rad)           ; cosine of latitude

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dxdy = dy * dx                  ; dydx(nlat)
   
   nlat2=dimsizes(lat2)
   nlon2=dimsizes(lon2)
   
   wgt  = new((/nlat2, nlon2/), typeof(ccl_sp({latS2:latN2},{lonL2:lonR2})))
   wgt  = conform(wgt, dxdy, 0)
    
   ccl_2 =wgt_areasum2(ccl_sp({latS2:latN2},{lonL2:lonR2}), wgt, 0)

   print(ccl_2)
   delete([/clat,dx,dy,dxdy,wgt/])

    ;3
   clat = cos(lat3 * rad)           ; cosine of latitude

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dxdy = dy * dx                  ; dydx(nlat)
   
   nlat3=dimsizes(lat3)
   nlon3=dimsizes(lon3)
   
   wgt  = new((/nlat3, nlon3/), typeof(ccl_sp({latS3:latN3},{lonL3:lonR3})))
   wgt  = conform(wgt, dxdy, 0)
    
   ccl_3 =wgt_areasum2(ccl_sp({latS3:latN3},{lonL3:lonR3}), wgt, 0)

   print(ccl_3)
   delete([/clat,dx,dy,dxdy,wgt/])

    ;4
   clat = cos(lat4 * rad)           ; cosine of latitude

   dx   = con * dlon * clat        ; dx at each latitude
   dy   = con * dlat               ; dy is constant
   dxdy = dy * dx                  ; dydx(nlat)
   
   nlat4=dimsizes(lat4)
   nlon4=dimsizes(lon4)
   
   wgt  = new((/nlat4, nlon4/), typeof(ccl_sp({latS4:latN4},{lonL4:lonR4})))
   wgt  = conform(wgt, dxdy, 0)
    
   ccl_4 =wgt_areasum2(ccl_sp({latS4:latN4},{lonL4:lonR4}), wgt, 0)

   ; print(ccl_4)
   ; delete([/clat,dx,dy,dxdy,wgt/])


   ; aa=(ccl_1+ccl_2+ccl_3+ccl_4)/ccl_global

   ; print(aa)



  



    ; ************************************************
    ; plot the correlations
    ; ************************************************
    wks  = gsn_open_wks("eps","global-ann-cummhw-space-90")             ; send graphics to PNG file
    plot = new (1, "graphic")
    setvalues NhlGetWorkspaceObjectId()
    "wsMaximumSize": 600000000
    ; 100000000
    end setvalues 
    


    ; gsn_define_colormap(wks,"cmocean_balance")     ;choose colormap
    res                       = True
    ; res@gsnMaximize           = True             ; make large
    res@gsnDraw               = False
    res@gsnFrame              = False
    res@cnFillOn              = True              ;turn on color
    res@cnLinesOn             = False             ;turn off contour lines
    res@cnLineLabelsOn        = False             ;turn off contour line labels

    ; res@mpLandFillColor         = "gray88"          ; make land white
    ; res@mpOutlineBoundarySets       = "AllBoundaries"
    ; res@mpDataBaseVersion      = "MediumRes" ;更好的地图效果
    res@mpGeophysicalLineThicknessF=0.001
    res@mpGeophysicalLineColor="white"
    ; For a grid this size, it is better to use raster mode. It will be 
    ; significantly faster, and will not go over NCL's 16mb default plot size.
    ; res@cnFillMode           = "RasterFill"       ; turn on raster mode

    res@vpHeightF          = 1.0               ; change aspect ratio of plot
    res@vpWidthF           = 0.7

    ; res@pmTickMarkDisplayMode = "Always"

    ; res@mpProjection      = "Robinson"        ; choose projection
    ; ; res@mpFillOn          = False
    ; res@mpGridAndLimbOn   = True       
    ; res@mpPerimOn         = False
    ; res@mpGridLatSpacingF =  90               ; change latitude  line spacing
    ; res@mpGridLonSpacingF = 180.              ; change longitude line spacing
    ; res@mpGridLineColor   = "transparent"     ; trick ncl into drawing
    ; res@vpHeightF          = 0.5               ; change aspect ratio of plot
    ; res@vpWidthF           = 0.6

    res@mpLandFillColor     = "white" 

    ; res@gsnAddCyclic          = True
    res@mpMinLatF              = -70                ; map area
    res@mpMaxLatF              = 70                 ;latitudes
    ; res@mpMinLonF              = lonL
    ; res@mpMaxLonF              = lonR
    res@mpCenterLonF          = 240.0
    ;---添加国界线
    ; res@mpOutlineBoundarySets = "National"

    ; res@mpFillOn              = True            ;turn off default background gray

    res@lbLabelBarOn          = False            ; turn off individual label bars
    ; res@lbOrientation = "vertical"

    res@cnLevelSelectionMode = "ManualLevels"      ; manually set cn levels

    res@cnLabelBarEndStyle   = "ExcludeOuterBoxes" ; remove the two end boxes from the labelbar
    res@pmLabelBarWidthF     = 0.63
    res@pmLabelBarHeightF     = 0.05
    res@pmLabelBarParallelPosF = 0.47  ;默认 0.5
    res@pmLabelBarOrthogonalPosF= 0.07
    res@tmXBLabelFontHeightF=0.01
    res@tmYLLabelFontHeightF=0.01
 
    res@tmXTOn                 = False    ; bottom off
    res@tmYROn                 = False    ; right off

    res@tmXBMajorOutwardLengthF = False
    res@tmXBMinorOutwardLengthF =False
    res@tmYLMajorOutwardLengthF = False
    res@tmYLMinorOutwardLengthF =False
    res@tmXBMajorLengthF=0.005
    res@tmXBMinorLengthF=0.00
    ; res@tmXTMajorLengthF
    res@tmYLMajorLengthF=0.005
    res@tmYLMinorLengthF=0.005
    ; res@gsnLeftString    ="(a) Frequency"
    ; ; res@gsnRightString    ="count/yr"
    ; res@gsnStringFontHeightF=0.023
    ; res@cnMinLevelValF       = 1.0                 ; min level
    ; res@cnMaxLevelValF       = 4.0                 ; max level
    ; res@cnLevelSpacingF      =  0.2                  ; contour level spacing
    ; res@cnFillPalette        = "MPL_YlOrRd"              ; choose colormap
    ; plot(0) = gsn_csm_contour_map(wks,ccl_sp(0,:,:),res)
    ; delete(res@cnMinLevelValF)
    ; delete(res@cnMaxLevelValF)
    ; delete(res@cnLevelSpacingF)
    ; res@tiMainString          = "(c) Cumulative Intensity"
    ; res@gsnCenterStringFontHeightF=0.015
    ; res@gsnCenterString    ="(c) Cumulative Intensity"
    ; res@gsnLeftString     = " "
    ; res@gsnRightString    = "~S~o~N~C x Days"
    ; res@tiMainOffsetYF        =  0.04                   ;-- move title upward
    ; res@gsnLeftString    ="Duration"
    ; res@gsnRightString    ="Days"
    res@cnMinLevelValF       = 0.0                 ; min level
    res@cnMaxLevelValF       = 150.0                 ; max level
    res@cnLevelSpacingF      =  20.0                  ; contour level spacing
    cmap=read_colormap_file("MPL_RdPu")
    res@cnFillPalette        = cmap(:,:)              ; choose colormap

    res@lbTitleOn        =  True                ; turn on title
  res@lbTitleString    = "(~S~o~N~C days)/counts"                ; title string
  res@lbTitlePosition  = "Right"              ; title position
  res@lbTitleFontHeightF= .011                ; make title smaller
  res@lbTitleDirection = "Across"             ; title direction

  res@lbTitleJust = "CenterRight"
  res@lbTitleExtentF = 0.001
  res@lbTitleOffsetF = -0.04
    plot(0) = gsn_csm_contour_map(wks,ccl_sp(:,:),res)
    delete(res@cnMinLevelValF)
    delete(res@cnMaxLevelValF)
    delete(res@cnLevelSpacingF)


    ; ;;添加子区域线条;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; resl               = True
    ; resl@gsLineColor   = "black"
    ; resl@gsLineThicknessF= 2.0
    ; resl@gsLineDashPattern=0

    ; ;atlnatic
    ; laty=(/5,45,45,5,5/)
    ; lonx=(/300,300,350,350,300/) 
    ; plot0 = gsn_add_polyline(wks,plot(0),lonx,laty,resl)

    ; ;tropical eastern pacific
    ; laty=(/-15,15,15,-15,-15/)
    ; lonx=(/210,210,280,280,210/) 
    ; plot1 = gsn_add_polyline(wks,plot(0),lonx,laty,resl)

  

    ; ;north pacific
    ; laty=(/20,60,60,20,20/)
    ; lonx=(/120,120,240,240,120/) 
    ; plot3 = gsn_add_polyline(wks,plot(0),lonx,laty,resl)

    

    

    ; ;Southwest Pacific
    ; laty=(/-58,-38,-38,-58,-58/)
    ; lonx=(/150,150,220,220,150/)  
    ; plot6 = gsn_add_polyline(wks,plot(0),lonx,laty,resl)

   



    draw(plot)

;   txres               = True                     ; text mods desired
;   txres@txFontHeightF = 0.01                     ; font smaller. default big

; ;                                      x ,y" etc.
; ;North Atlantic
;   gsn_text_ndc(wks,"North Atlantic",0.69,0.495,txres) 

;   ;Southwest Pacific
;   gsn_text_ndc(wks,"Southwest Pacific",0.40,0.375,txres) 

; ;NINO Pacific
;   gsn_text_ndc(wks,"Tropical Eastern Pacific",0.52,0.46,txres) 
; ;Westnorth Pacific
;   gsn_text_ndc(wks,"North Pacific ",0.36,0.53,txres) 





    frame(wks)

