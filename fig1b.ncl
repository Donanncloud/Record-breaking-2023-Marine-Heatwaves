    latS   =  -60
    latN   =60
    lonL   =  0
    lonR   =  360

    ;================================================
    diri  = "/public/home/dongty/code/2023global/"
    diri2="/public/home/dongty/code/2023global/ann/figureV2/"
    diri3  = "/public/home/dongty/code/2023global/ann/cate/"
    
    fili1 = "OISST_4mhw_counts_durs_mean-cum_global_allyear_90th_1982to2023.nc"              
    fili2 = "cate_fraction-time-90.nc"        
    fili3 = "OISST_year_category_cates_global.nc"     
   


    a1    = addfile(diri+fili1,"r")
    a2    = addfile(diri2+fili2,"r")
    a3    = addfile(diri3+fili3,"r")
    ; a4    = addfile(diri+fili4,"r")

    ccl     = a1->ccl(:,:,{latS:latN},{lonL:lonR})
    cate_fract      = (a2->cate_fract)
    counts     = a3->counts(:,{latS:latN},{lonL:lonR})

  counts:=where(counts.eq."nan",-9.96921e+36,counts)
  counts@_FillValue=-9.96921e+36
  counts:=where(counts.ge.0.0 , 1.0,counts@_FillValue)

     lat=a1->lat({latS:latN})
    lon=a1->lon({lonL:lonR})


    nlat=dimsizes(lat)
    nlon=dimsizes(lon)
    lat@units="degrees-north"
    lon@units="degrees-east"

  ; durs:=where(durs.eq."nan",-9.96921e+36,durs)
  ; durs@_FillValue=-9.96921e+36

  ; printVarSummary(durs)

  ; counts=dim_sum_n( durs(:,:,5:7,:,:),2)
  ; counts:=dim_avg_n(counts,1)
  ; counts=where(counts.eq.0, -9.96921e+36 ,counts)
  counts!1="lat"
  counts!2="lon"
  counts&lat=lat
  counts&lon=lon
  counts&lat@units="degrees-north"
  counts&lon@units="degrees-east"
 

    ; print(cate_fract)
    ; inten_means     = a3->inten_means
    ; inten_cums     = a4->inten_cums



    year=ispan(1982,2023,1)
    nyear=42



    ccl!1="year"
    ccl!2="lat"
    ccl!3="lon"
    ccl&year=year
    ccl&lat=lat
    ccl&lon=lon
    ccl&lat@units="degrees-north"
    ccl&lon@units="degrees-east"
    printVarSummary(ccl)


    ; ;space 
    ; ccl_sp=dim_avg_n(ccl, 1)

    ; ccl_sp!1="lat"
    ; ccl_sp!2="lon"
    ; ccl_sp&lat=lat
    ; ccl_sp&lon=lon
    ; ccl_sp&lat@units="degrees-north"
    ; ccl_sp&lon@units="degrees-east"

    allcum=ccl(0,:,:,:)*ccl(3,:,:,:)

    allcum!0="year"
    allcum!1="lat"
    allcum!2="lon"
    allcum&year=year
    allcum&lat=lat
    allcum&lon=lon
    allcum&lat@units="degrees-north"
    allcum&lon@units="degrees-east"
    printVarSummary(allcum)

    area=allcum
    area:=where(allcum.gt.0.0,1.0,allcum)

    ; print(allcum(1:10,100,100))
    ; print(area(1:10,100,100))


    ; time 区域平均
    rad    = 3.14159/180.0
    jlat   = dimsizes(lat)
    dlon   = abs(lon(2)-lon(1))
    dx     = dlon*cos(lat*rad)
    dy     = new (jlat, typeof(dx))
    dy(0)  = abs(lat(2)-lat(1))
    dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*0.5   
    dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))
    dxdy  = dx*dy
    copy_VarCoords(allcum(0,:,:),dxdy)

    printVarSummary(dxdy)
    copy_VarCoords(ccl(0,0,:,:),dxdy)

    allcum_ts =wgt_areaave_Wrap(allcum, dxdy, 1.0, 1)

    meanint_ts=wgt_areaave_Wrap(ccl(2,:,:,:), dxdy, 1.0, 1)
    durs_ts=wgt_areaave_Wrap(ccl(0,:,:,:)*ccl(1,:,:,:), dxdy, 1.0, 1)

  
    ; ;    乘上发生热浪的覆盖范围
    ; allcum_ts:=allcum_ts(:)*cate_fract
    ; print(cate_fract)
    ; print(meanint_ts)

    delete([/rad,jlat,dlon,dx,dy,dxdy/])

    re   = 6.37122e03
    rad  = 4.0 * atan(1.0) / 180.0
    con  = re * rad                 
    clat = cos(lat * rad)           ; cosine of latitude

    dlon = (lon(2) - lon(1))        ; assume dlon is constant
    dlat = (lat(2) - lat(1))        ; assume dlat is constant

    dx   = con * dlon * clat        ; dx at each latitude
    dy   = con * dlat               ; dy is constant
    dxdy = dy * dx                  ; dydx(nlat)




    wgt  = new((/nlat, nlon/), typeof(allcum))
    wgt  = conform(wgt, dxdy, 0)
   
    counts_all= wgt_areasum2(counts(:,:,:), wgt, 0)/10^8

    areaall= (wgt_areasum2(area(:,:,:), wgt, 0)/10^8 ) /counts_all *100
   
    
    ; print(counts_all)
    ; print(meanint_ts)
    ; print(avg(durs_ts(:40)))
    ; print(durs_ts)

    ; print(areaall)

    ; allcum_ts!0="year"
    ; allcum_ts&year=year
    ; printVarSummary(allcum_ts) 
    
    ;regress 
    rc_durs = regline_stats(year, durs_ts)                ; degC/year
    rc_durs@long_name = "trend"
    rc_durs_plot= rc_durs@Yest                           ; regression values 


    rc_mean_int = regline_stats(year, meanint_ts)                ; degC/year
    rc_mean_int@long_name = "trend"
    rc_mean_plot= rc_mean_int@Yest                           ; regression values 

    rc_area= regline_stats(year, areaall)                ; degC/year
    rc_area@long_name = "trend"
    rc_area_plot= rc_area@Yest                           ; regression values 


   ;average
   value_durs=rc_durs@Yavg
   value_mean_int=rc_mean_int@Yavg
   value_area=rc_area@Yavg
 
 rc=regline(year, areaall)
 df   = rc@nptxy-2
 prob =  betainc(df/(df+rc@tval^2), df/2.0, 0.5) 
print(rc)
print(prob)
; print(rc_area@pval)
  


    ;----------------------------------------------------------------------------------------
    ;画图
    ;----------------------------------------------------------------------------------------
    wks   = gsn_open_wks ("pdf", "3lines-v4-fig1b" )  ; send graphics to PNG file

    ; left variable
    resL           = True 

    resL@gsnFrame  = False           ; Don't draw plot or advance frame
    resL@gsnDraw   = False

    ; resL@vpXF      = 0.15
    resL@vpWidthF  = 0.55            ; Set width and height of plot.
    resL@vpHeightF = 0.3

    resL@trXMinF = 1980.                 ; Set min/max of X axis.
    resL@trXMaxF = 2025.
    resL@trYMinF                = 0
    resL@trYMaxF                = 125

    xindexc = (/1980,1985,1990,1995,2000,2005,2010,2015,2020,2023,2025/)
    xlabelc = (/"1980"," ","1990"," ","2000"," ","2010","2015"," ","2023"," "/)
    resL@tmXBMode             = "Explicit"
    ;以下两行代码都要以上面一行代码为前提，否则………………
    resL@tmXBValues            = xindexc  
    resL@tmXBLabels            = xlabelc

    resL@xyLineColor       = "black"   ; To distinguish the lines
    resL@xyLineThicknesses = 4.0      ; thicker line

    ; tickmarks:
    resL@tmXTOn              = False   ; Turn off top tickmarks

    resL@tmYLMinorOn         = False   ; Turn off left minor tickmarks

    resL@tmXBMode            = "Explicit"
    resL@tmXBLabelFontHeightF = 0.01



    resL@tmXTOn                 = False    ; bottom off
    resL@tmXBMajorOutwardLengthF = False
    resL@tmXBMinorOutwardLengthF =False
    resL@tmYLMajorOutwardLengthF = False
    resL@tmYLMinorOutwardLengthF =False
    resL@tmYRMajorOutwardLengthF = False
    resL@tmYRMinorOutwardLengthF =False
    resL@tmXBMajorLengthF=0.01

    ; res@tmXTMajorLengthF
    resL@tmYLMajorLengthF=0.01
    resL@tmYLMinorLengthF=0.005
    resL@tmYRMajorLengthF=0.01
    resL@tmYRMinorLengthF=0.005

    ; resL@tiMainString                = "Four curves with four separate Y axes"

    resL@tiYAxisString               = "Duration (days)"  ; axis string
    resL@tiYAxisFontColor            = resL@xyLineColor
    resL@tiYAxisFontHeightF          = 0.015
    resL@tmYLLabelFontColor=resL@xyLineColor
    resL@xyDashPatterns         = 0 


    ; right variable
    resR                        = True                      
    resR@xyLineColor            = "darkgoldenrod2"
    resR@xyDashPatterns         = 0                   ; dashed line for 2nd
    resR@xyLineThicknesses      = 4.0                   ; thicker line

    resR@trXMinF = resL@trXMinF
    resR@trXMaxF = resL@trXMaxF
    resR@trYMinF                = 0.5
    resR@trYMaxF                = 1.5

    resR@tmYRLabelFontColor=resR@xyLineColor

    resR@tiYAxisString    = "Mean intensity (~S~o~N~C)"   ; axis string
    resR@tiYAxisFontColor = resR@xyLineColor
    resR@tiYAxisOffsetXF=-0.005


    resR@tmYRMinorOn         = False   ; Turn off right minor tickmarks

    ; Laplacian variable
    resRa                        = True

    resRa@trXMinF                = resL@trXMinF
    resRa@trXMaxF                = resL@trXMaxF
    
    resRa@trYMinF                = 40
    resRa@trYMaxF                = 100

    resRa@xyDashPatterns         = 0                   ; dashed line for 2nd
    resRa@xyLineThicknesses      = 4.0                   ; thicker line
    resRa@xyLineColor            = "deepskyblue3"

    resRa@tiYAxisString    = "Percentage of Area (%)"   ; axis string
    resRa@tiYAxisFontColor = resRa@xyLineColor
    resRa@tmYRLabelFontColor=resRa@xyLineColor

    resRa@tmYRMinorOn = False   ; Turn off right minor tickmarks

    resRa@amOrthogonalPosF = 0.65 ; Move "l" axis to right a little

    resL@tmXBLabelFontHeightF = 0.012
    resL@tmYLLabelFontHeightF = 0.012
    resR@tmYRLabelFontHeightF = 0.012
    resRa@tmYRLabelFontHeightF = 0.012

    ; resL@tmBorderThicknessF=0.2
    resR@tmBorderThicknessF=0.0
    resRa@tmBorderThicknessF=0.0
    
    resL@tmYLMinorOn=True
    resR@tmYRMinorOn=True
    resRa@tmYRMinorOn=True
    resL@tmYMajorGrid                = False
  ; resL@tmXMajorGridThicknessF      = 1.0
    resL@tmYMajorGridThicknessF      = 1.0
  ; resL@tmXMajorGridLineDashPattern = 2
    resL@tmYMajorGridLineDashPattern = 2
    ; Create the plot with the three variables and three different axes.
    ; Nothing is drawn at this point.
    ;Left


    resL@gsnYRefLineThicknesses = 2.
    resL@gsnYRefLineDashPattern=2
    resL@gsnYRefLineColors   =(/"black"/)
    resL@gsnYRefLine           = sprintf("%6.2f", value_durs)             ; create a reference line 

    resR@gsnYRefLineThicknesses = 2.
    resR@gsnYRefLineDashPattern=2
    resR@gsnYRefLineColors   =(/"darkgoldenrod2"/)
    resR@gsnYRefLine           = sprintf("%6.2f", value_mean_int)             ; create a reference line 

    ; resRa@gsnYRefLineThicknesses = 3.
    ; resRa@gsnYRefLineDashPattern=12
    ; resRa@gsnYRefLineColors   =(/"deepskyblue2"/)
    ; resRa@gsnYRefLine           = sprintf("%6.2f", value_area)             ; create a reference line  
    ; print(value_area)
    ; print(value_mean_int)
    ; print(value_durs)
    plot = gsn_csm_xy3(wks,year,durs_ts(:),meanint_ts(:),areaall(:),resL,resR,resRa)
    
    ; plot2=gsn_csm_xy3(wks,year,durs_ts(:),meanint_ts(:),areaall(:),resL,resR,resRa)


;-------------------------------------
; 手动添加蓝色 RefLine for resRa
;-------------------------------------
x_min = 1980.0
x_max = 2025.0
y_val = value_area

xx = (/ x_min, x_max /)
yy = (/ y_val, y_val /)

res_line = True
res_line@gsLineColor = "deepskyblue2"
res_line@gsLineThicknessF = 2.0
res_line@gsLineDashPattern = 2

dum_line = gsn_add_polyline(wks, plot, xx, yy, res_line)

  ; res_lines                   = True                  
  ; res_lines@gsLineDashPattern = 0                  
  ; res_lines@gsLineThicknessF  = 3.0                  
 
  res_line@gsLineColor       = "grey36"                 
  xx = (/2006,2009/)
  yy = (/117,117/)

  dum_1a = gsn_add_polyline(wks, plot,xx,yy,res_line) 

;;add trend
    
  res_text               = True
  res_text@txFontHeightF = 0.012
  
  res_text@txJust             = "CenterLeft"  

                  
  res_text@txFontColor   = "grey36"   
  dum_1= gsn_add_text(wks, plot,"Trend:",1981,117,res_text)
  dum_1b= gsn_add_text(wks, plot,"Average:",2000,117,res_text)
 

  res_text@txFontColor   = "black"            
  dum_1= gsn_add_text(wks, plot,sprintf("%6.2f", rc_durs*10)+"*"+" days"+" decade"+"~S~-1~N~",1985,117,res_text)
  dum_11= gsn_add_text(wks, plot,sprintf("%4.2f", value_durs)+" days",2017,40,res_text)

  res_text@txFontColor   = "darkgoldenrod2"  
  dum_2= gsn_add_text(wks, plot,sprintf("%6.2f", rc_mean_int*10)+"*"+" ~S~o~N~C"+" decade"+"~S~-1~N~",1985,107,res_text)
  dum_22= gsn_add_text(wks, plot,sprintf("%4.2f", value_mean_int)+" ~S~o~N~C",1982,64,res_text)
  
  res_text@txFontColor   = "deepskyblue3"  
  dum_3= gsn_add_text(wks, plot,sprintf("%6.2f", rc_area*10)+"*"+" %"+" decade"+"~S~-1~N~",1985,97,res_text)
  dum_33= gsn_add_text(wks, plot,sprintf("%4.2f", value_area)+" %",1982,78,res_text)





     draw(plot)
     frame(wks)