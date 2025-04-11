     latS   =  -60
     latN   =60
     lonL   =  0
     lonR   =  360

     ;================================================

     diri2="/public/home/dongty/code/2023global/"
     fili1 = "OISST_4mhw_counts_durs_mean-cum_global_allyear_90th_1982to2023.nc"              


     a1    = addfile(diri2+fili1,"r")


     cum0     = a1->ccl(3,:,{latS:latN},{lonL:lonR})
     fre     = a1->ccl(0,:,{latS:latN},{lonL:lonR})

     allcum=fre*cum0
     ; cate_fract      = a2->cate_fract
     ; inten_means     = a3->inten_means
     ; inten_cums     = a4->inten_cums

     lat=a1->lat({latS:latN})
     lon=a1->lon({lonL:lonR})


     nlat=dimsizes(lat)
     nlon=dimsizes(lon)
     lat@units="degrees-north"
     lon@units="degrees-east"

     year=ispan(1982,2023,1)
     nyear=42



     allcum!0="year"
     allcum!1="lat"
     allcum!2="lon"
     allcum&year=year
     allcum&lat=lat
     allcum&lon=lon
     allcum&lat@units="degrees-north"
     allcum&lon@units="degrees-east"
     printVarSummary(allcum)



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
   

   
   wgt  = new((/nlat, nlon/), typeof(allcum))
   wgt  = conform(wgt, dxdy, 0)
    
    allcum_ts =wgt_areasum2(allcum(:,:,:), wgt, 0)
    allcum_ts:=allcum_ts/10^10

    allcum_ts!0="year"
    allcum_ts&year=year
    printVarSummary(allcum_ts) 

    print(allcum_ts)



     cum=allcum_ts
     print(cum)
     mean=dim_avg(cum({1982:2022}))
     std_1=dim_stddev_n_Wrap(cum(:), 0) 
     print(mean)
     print(std_1)

     cum_2023=cum(40:41)
     year3=ispan(2022,2023,1)


     mean_std_1     = new ((/2,46/), float)
     ; mean_std_1(0,:) = mean-std_1
     mean_std_1(0,:) = mean+std_1
     mean_std_1(1,:) = mean
     
     year2=ispan(1980,2025,1)
     mean_std_1!1="year2"
     mean_std_1&year2=year2
     printVarSummary(mean_std_1)

     mean_std_2     = new ((/2,46/), float)
     ; mean_std_2(0,:) = mean-std_1*2
     mean_std_2(0,:) = mean+std_1*2
     mean_std_2(1,:) = mean

     copy_VarMeta(mean_std_1, mean_std_2)

     mean_std_3     = new ((/2,46/), float)
     ; mean_std_3(0,:) = mean-std_1*3
     mean_std_3(0,:) = mean+std_1*3
     mean_std_3(1,:) = mean

     copy_VarMeta(mean_std_1, mean_std_3)

     aa=cum(41)/mean
     print(aa)

     ; ************************************************
     ; plot the correlations
     ; ************************************************
     wks  = gsn_open_wks("eps","cum*area*fre-ann-global-time-90-fig1a")             ; send graphics to PNG file
     plot = new (2, "graphic")
     plot_shad = new (3, "graphic")

     res               = True                   ; plot mods desivioletred
     res@gsnDraw            = False             ; don't draw yet
     res@gsnFrame           = False             ; don't advance frame yet
     
     ;res@gsnPaperOrientation      ="portrait"
     res@gsnMaximize               = True

     res@vpHeightF          = 0.4               ; change aspect ratio of plot
     res@vpWidthF           = 0.7

     res@trXMinF=1980
     res@trXMaxF=2025
     res@trYMinF=0.0
     res@trYMaxF=6.0  
     ;=====================================================================
     xindexc = (/1980,1985,1990,1995,2000,2005,2010,2015,2020,2023,2025/)
     xlabelc = (/"1980"," ","1990"," ","2000"," ","2010","2015"," ","2023"," "/)

     ;这里设置x坐标轴的具体表示，而非缺省值
     res@tmXBMode             = "Explicit"
     ;以下两行代码都要以上面一行代码为前提，否则………………
     res@tmXBValues            = xindexc  
     res@tmXBLabels            = xlabelc
     res@tmXBLabelFontHeightF = 0.018
     res@tmYLLabelFontHeightF = 0.018
     ; ------------------------------------------------------

     res@xyMonoDashPattern      = True     ; force all solid lines
     ;---Turn off some tickmarks
     res@tmXTOn                 = False    ; bottom off
     res@tmYROn                 = False    ; right off

     res@tmXBMajorOutwardLengthF = False
     res@tmXBMinorOutwardLengthF =False
     res@tmYLMajorOutwardLengthF = False
     res@tmYLMinorOutwardLengthF =False
     res@tmXBMajorLengthF=0.01

     ; res@tmXTMajorLengthF
     res@tmYLMajorLengthF=0.01
     res@tmYLMinorLengthF=0.005
     ; res@tmYRMajorLengthF

     res@gsnRightString ="(10~S~10~N~ ~S~o~N~C days km~S~2~N~)"
     res@gsnRightStringFontHeightF=0.02
     res@gsnRightStringFont   = "helvetica"

     res@tiYAxisFontHeightF=0.02
     res@tiXAxisFontHeightF=0.02
     res@tiYAxisString  ="MHW Activity"
     ; res@tiXAxisString  ="Year"
     res@tmYRLabelsOn         = False           ; turn on right axis labels
     ;设置图例
     ; res@pmLegendDisplayMode    = "Always"       ;打开图例     ; turn on legend
     ; res@pmLegendWidthF         = 0.0                ; Change width and
     ; res@pmLegendHeightF        = 0.18                ; height of legend.
     ; res@lgPerimOn              = False   ; 这个是在图例周围加一圈框框
     ; res@lgLabelFontHeightF     = .02    ; 标签字体高度

     res@xyLineThicknesses      = 9.0 
     ;以下是画三个图，中间注意要把前一张设置好的参数delete掉，这样能释放内存，也不容易出错

     ; res@tiMainFontHeightF  =0.03      ;设置标题字体的高度 
     ; res@gsnRightStringFontHeightF  =0.027      ;设置标题字体的高度 
     res@xyLineColors           = (/"deepskyblue3"/)  ;使用colors1来画第一张图的颜色
     res@xyMarkerColors= (/"deepskyblue3"/) 
     res@xyMarkLineModes = (/"MarkLines"/)
     res@xyMarkers           = (/16/)
     ; res@gsnLeftString    ="(e) Frequency"
     ; res@gsnRightString    ="count/yr"

     ; res@trXMinF=1982
     ; res@trXMaxF=2024

     ; plot(0)=gsn_csm_xy(wks,year,ccl_ts(0,:),res)

     res@gsnYRefLine           = mean             ; create a reference line  
     res@gsnYRefLineDashPattern=12
     res@gsnYRefLineThicknesses = 5.
     res@gsnYRefLineColors   =(/"deepskyblue2"/)
     ; res@gsnLeftString    ="(a) Cumulative Intensity"
     ; res@gsnRightString    ="~S~o~N~C x Days" 
     res@tmYLMinorOn         =True
     plot(0)=gsn_csm_xy(wks,year(0:40),cum(0:40),res)

     res@xyLineColors           = (/"violetred"/)
     res@xyMarkerColors= (/"violetred"/)
     delete(res@gsnRightString)
     plot(1)=gsn_csm_xy(wks,year3,cum_2023(:),res)



     delete([/res@xyMarkLineModes,res@xyLineColors/])



     res@gsnXYFillColors = "lightskyblue"
     res@gsnXYFillOpacities=0.4
     res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
     plot_shad(0)  = gsn_csm_xy (wks,year2,mean_std_1(0:1,:),res)  ; Create filled XY plot.

     res@gsnXYFillOpacities=0.2
     ; res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
     plot_shad(1)  = gsn_csm_xy (wks,year2,mean_std_2(0:1,:),res)  ; Create filled XY plot.

    
     res@gsnXYFillOpacities=0.1
     ; res@xyLineColor     = -1                           ; We don't want the line, so make it transparent.
     plot_shad(2)  = gsn_csm_xy (wks,year2,mean_std_3(0:1,:),res)  ; Create filled XY plot.
     
     overlay(plot(1),plot(0))


     overlay(plot(0),plot_shad(0))
     overlay(plot(0),plot_shad(1)) 
     overlay(plot(0),plot_shad(2)) 

     ;=================================================
     ; add greek symbol to plot with gsn_text
     ;=================================================
     ; txres               = True                      ; text mods desivioletred
     ; txres@txFontHeightF = 0.07                      ; text font height

     ; gsn_text_ndc(wks,"xxxx",50,5,txres)

     ;  gsres               = True
     ;  gsres@gsMarkerIndex = 14          ; Use filled dots for markers.

     ;  gsres@gsMarkerColor      = "violetred"
     ;  gsres@gsMarkerThicknessF = 28.0
     ; ; gsn_polymarker(wks,plot,0.95,0.645,gsres)
     ;  gsn_polymarker_ndc(wks, 0.955,0.687,gsres) ;x,y

     txres               = True                     ; text mods desivioletred
     txres@txFontHeightF = 0.02                     ; font smaller. default big
     txres@txFontColor  ="dodgerblue4"
     ; txres@gsnXYFillOpacities=0.9
     ; txres@cnLineDrawOrder=True
     gsn_text_ndc(wks,"+~F21~1~F33~s",0.2,0.45,txres) 
     txres@txFontColor  ="dodgerblue2"
     gsn_text_ndc(wks,"+~F21~2~F33~s",0.2,0.55,txres) 
     txres@txFontColor  ="deepskyblue3"
     ; txres@gsnXYFillOpacities=0.1
     gsn_text_ndc(wks,"+~F21~3~F33~s",0.2,0.65,txres) 

     gsn_text_ndc(wks,"1982-2022 mean",0.85,0.38,txres) 
     txres@txFontColor  ="violetred"
     txres@txFontHeightF = 0.03 
     gsn_text_ndc(wks,"2023",0.88,0.73,txres) 
     txres@txFontHeightF = 0.025
     gsn_text_ndc(wks,"(5.36)",0.88,0.69,txres)  
     
     txres@txFontHeightF = 0.025
     gsn_text_ndc(wks,"~F37~Z~",0.941,0.73,txres)  



     draw(plot)
     frame(wks)


