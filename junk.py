#print(tab_pe_SAMI.columns)
for cols in (tab_pe_SAMI.columns):
    tab_pe_SAMI[cols].units=''
    tab_pe_SAMI[cols].description=''
#print(tab_pe_SAMI[tab_pe_SAMI['name']=='manga-11757-9102'])

tab_pe_SAMI_select=tab_pe_SAMI['name','lSFR','z_gas','z_stars','FoV','nx','ny','Re_kpc','log_Mass_corr','error_Mass','e_lSFR','log_age_mean_LW','s_log_age_mean_LW','log_ZH_mean_LW','s_log_ZH_mean_LW','Av_ssp_stats_mean','Av_ssp_stats_stddev','Av_w1','e_Av_w','log_Mass','log_SFR_ssp','log_NII_Ha_cen_mean','log_NII_Ha_cen_stddev','log_OIII_Hb_cen_mean','log_OIII_Hb_cen_stddev','log_SII_Ha_cen_mean','log_SII_Ha_cen_stddev','log_OII_Hb_cen_mean','log_OII_Hb_cen_stddev','EW_Ha_cen_mean','EW_Ha_cen_stddev','ZH_LW_Re_fit','e_ZH_LW_Re_fit','alpha_ZH_LW_Re_fit','e_alpha_ZH_LW_Re_fit','ZH_MW_Re_fit','e_ZH_MW_Re_fit','alpha_ZH_MW_Re_fit','e_alpha_ZH_MW_Re_fit','Age_LW_Re_fit','e_Age_LW_Re_fit','alpha_Age_LW_Re_fit','e_alpha_Age_LW_Re_fit','Age_MW_Re_fit','e_Age_MW_Re_fit','alpha_Age_MW_Re_fit','e_alpha_Age_MW_Re_fit','Re_arc','DL','DA','P.A.','Ellipticity','Inclination','elip_ab','log_Mass_gas','rat_vel_sigma','e_rat_vel_sigma','log_SFR_SF','log_SFR_D_C','OH_O3N2_cen','e_OH_O3N2_cen','OH_N2_cen','e_OH_N2_cen','OH_ONS_cen','e_OH_ONS_cen','OH_R23_cen','e_OH_R23_cen','OH_pyqz_cen','e_OH_pyqz_cen','OH_t2_cen','e_OH_t2_cen','OH_M08_cen','e_OH_M08_cen','OH_T04_cen','e_OH_T04_cen','OH_dop_cen','e_OH_dop_cen','OH_O3N2_EPM09_cen','e_OH_O3N2_EPM09_cen','log_OI_Ha_cen','e_log_OI_Ha_cen','Ha_Hb_cen','e_Ha_Hb_cen','log_NII_Ha_Re','e_log_NII_Ha_Re','log_OIII_Hb_Re','e_log_OIII_Hb_Re','log_SII_Ha_Re','e_log_SII_Ha_Re','log_OII_Hb_Re','e_log_OII_Hb_Re','log_OI_Ha_Re','e_log_OI_Ha_Re','EW_Ha_Re','e_EW_Ha_Re','Ha_Hb_Re','e_Ha_Hb_Re','log_NII_Ha_ALL','e_log_NII_Ha_ALL','log_OIII_Hb_ALL','e_log_OIII_Hb_ALL','log_SII_Ha_ALL','e_log_SII_Ha_ALL','log_OII_Hb_ALL','e_log_OII_Hb_ALL','log_OI_Ha_ALL','e_log_OI_Ha_ALL','EW_Ha_ALL','e_EW_Ha_ALL','Ha_Hb_ALL','Sigma_Mass_cen','e_Sigma_Mass_cen','Sigma_Mass_Re','e_Sigma_Mass_Re','Sigma_Mass_ALL','e_Sigma_Mass_ALL','T30','ZH_T30','ZH_Re_T30','a_ZH_T30','T40','ZH_T40','ZH_Re_T40','a_ZH_T40','T50','ZH_T50','ZH_Re_T50','a_ZH_T50','T60','ZH_T60','ZH_Re_T60','a_ZH_T60','T70','ZH_T70','ZH_Re_T70','a_ZH_T70','T80','ZH_T80','ZH_Re_T80','a_ZH_T80','T90','ZH_T90','ZH_Re_T90','a_ZH_T90','T95','ZH_T95','ZH_Re_T95','a_ZH_T95','T99','ZH_T99','ZH_Re_T99','a_ZH_T99','log_Mass_gas_Av_gas_OH','log_Mass_gas_Av_ssp_OH','vel_ssp_2','e_vel_ssp_2','vel_Ha_2','e_vel_Ha_2','vel_ssp_1','e_vel_ssp_1','vel_Ha_1','e_vel_Ha_1','log_SFR_ssp_100Myr','log_SFR_ssp_10Myr','vel_disp_Ha_cen','vel_disp_ssp_cen','vel_disp_Ha_1Re','vel_disp_ssp_1Re','vel_disp_Ha_1Re_mean','vel_disp_ssp_1Re_mean','KIN_Ha_1','e_KIN_Ha_1','KIN_ssp_1','e_KIN_ssp_1','KIN_Ha_05','e_KIN_Ha_05','KIN_ssp_05','e_KIN_ssp_05','log_Mass_in_Re','ML_int','ML_avg','F_Ha_cen','e_F_Ha_cen','Re_kpc_Mass','R50_kpc_V','e_R50_kpc_V','R50_kpc_Mass','e_R50_kpc_Mass','log_Mass_corr_in_R50_V','e_log_Mass_corr_in_R50_V','log_Mass_gas_Av_gas_log_log']

print(tab_pe_SAMI_select[tab_pe_SAMI_select['name']=='manga-11834-9102'])
print(tab_pe_SAMI_select[tab_pe_SAMI_select['name']=='manga-11941-1902'])
#print(tab_all_select[tab_all_select['name']=='manga-11834-9102'])
#print(tab_pe_SAMI_select[tab_pe_SAMI_select['name']=='manga-11757-9102'])

tab_pe_SAMI_select.rename_column('lSFR','log_SFR_Ha')
tab_pe_SAMI_select.rename_column('e_lSFR','e_log_SFR_Ha')
tab_pe_SAMI_select.rename_column('log_NII_Ha_cen_mean','log_NII_Ha_cen')
tab_pe_SAMI_select.rename_column('log_NII_Ha_cen_stddev','e_log_NII_Ha_cen')  
tab_pe_SAMI_select.rename_column('log_OIII_Hb_cen_mean','log_OIII_Hb_cen')
tab_pe_SAMI_select.rename_column('log_OIII_Hb_cen_stddev','e_log_OIII_Hb_cen')
tab_pe_SAMI_select.rename_column('log_SII_Ha_cen_mean','log_SII_Ha_cen')
tab_pe_SAMI_select.rename_column('log_SII_Ha_cen_stddev','e_log_SII_Ha_cen')
tab_pe_SAMI_select.rename_column('log_OII_Hb_cen_mean','log_OII_Hb_cen')
tab_pe_SAMI_select.rename_column('log_OII_Hb_cen_stddev','e_log_OII_Hb_cen')
tab_pe_SAMI_select.rename_column('EW_Ha_cen_mean','EW_Ha_cen')
tab_pe_SAMI_select.rename_column('EW_Ha_cen_stddev','e_EW_Ha_cen')
tab_pe_SAMI_select.rename_column('rat_vel_sigma','vel_sigma_Re')# float64  
tab_pe_SAMI_select.rename_column('e_rat_vel_sigma','e_vel_sigma_Re')# float64  

tab_pe_SAMI_select.remove_column('log_age_mean_LW')
tab_pe_SAMI_select.remove_column('s_log_age_mean_LW')
tab_pe_SAMI_select.remove_column('log_ZH_mean_LW')
tab_pe_SAMI_select.remove_column('s_log_ZH_mean_LW')
tab_pe_SAMI_select.remove_column('Av_ssp_stats_mean')
tab_pe_SAMI_select.remove_column('Av_ssp_stats_stddev')
tab_pe_SAMI_select.remove_column('Av_w1') 
tab_pe_SAMI_select.remove_column('e_Av_w')
tab_pe_SAMI_select.remove_column('elip_ab')
tab_pe_SAMI_select.remove_column('Inclination')
tab_pe_SAMI_select.remove_column('z_gas')# float64  
tab_pe_SAMI_select.remove_column('z_stars')# float64  
tab_pe_SAMI_select.remove_column('nx')# int64  
tab_pe_SAMI_select.remove_column('ny')# int64  
tab_pe_SAMI_select.remove_column('log_Mass_corr')# float64  
tab_pe_SAMI_select.remove_column('vel_disp_Ha_1Re_mean') #float64  
tab_pe_SAMI_select.remove_column('vel_disp_ssp_1Re_mean')# float64  
tab_pe_SAMI_select.remove_column('KIN_Ha_1')# float64  
tab_pe_SAMI_select.remove_column('e_KIN_Ha_1')# float64  
tab_pe_SAMI_select.remove_column('KIN_ssp_1')# float64  
tab_pe_SAMI_select.remove_column('e_KIN_ssp_1')# float64  
tab_pe_SAMI_select.remove_column('KIN_Ha_05')# float64  
tab_pe_SAMI_select.remove_column('e_KIN_Ha_05')# float64  
tab_pe_SAMI_select.remove_column('KIN_ssp_05')# float64  
tab_pe_SAMI_select.remove_column('e_KIN_ssp_05')# float64 
tab_pe_SAMI_select.remove_column('Re_kpc_Mass')


tab_pe_SAMI_select.rename_column('error_Mass','e_log_Mass')# float64  

tab_val_select=tab_val_Re['name','Av_gas','e_Av_gas','Av_ssp','e_Av_ssp','Lambda','e_Lambda']
tab_pe_SAMI_select=join(tab_pe_SAMI_select,tab_val_select,keys=['name'],join_type='left')
tab_pe_SAMI_select.rename_column('Lambda','Lambda_Re')
tab_pe_SAMI_select.rename_column('e_Lambda','e_Lambda_Re')
tab_pe_SAMI_select.rename_column('Av_gas','Av_gas_Re')
tab_pe_SAMI_select.rename_column('e_Av_gas','e_Av_gas_Re')
tab_pe_SAMI_select.rename_column('Av_ssp','Av_ssp_Re')
tab_pe_SAMI_select.rename_column('e_Av_ssp','e_Av_ssp_Re')
tab_pe_SAMI_select.rename_column('P.A.','PA')
tab_pe_SAMI_select.rename_column('Ellipticity','ellip')                                                 
for cols in (tab_pe_SAMI_select.columns):
    tab_pe_SAMI_select[cols].units=''
    tab_pe_SAMI_select[cols].description=''

tab_pe_SAMI_info = ascii.read(dir_DR17+"tables/get_proc_elines_MaNGA.select.txt",delimiter=',', guess=True)#, fill_values=" ")  
for cols in (tab_pe_SAMI_select.columns):
    mask = (tab_pe_SAMI_info['NAME']==cols)
    tab_tmp=tab_pe_SAMI_info['DESC'][mask]
    if(len(list(tab_tmp))>0):
        desc_now=list(tab_tmp)[0]
        tab_tmp=tab_pe_SAMI_info['UNITS'][mask]
        units_now=list(tab_tmp)[0]
        #print(cols,',',units_now,',',desc_now)
        tab_pe_SAMI_select[cols].units=units_now
        tab_pe_SAMI_select[cols].description=desc_now

        
for cols in (tab_pe_SAMI_select.columns):
    #
    # Oxygen abundances
    #
    if (cols.find("OH")>-1):
        tab_pe_SAMI_select[cols].units='dex'
        if (cols.find("alpha")>-1):
            if (cols.find("e_OH")>-1):
                cal=cols.replace('e_OH_','')
                cal=cal.replace('_alpha_fit','')            
                tab_pe_SAMI_select[cols].description='Error in the slope of the O/H gradient using the calibrator '+cal
            else:
                cal=cols.replace('OH_','')
                cal=cal.replace('_alpha_fit','')
                tab_pe_SAMI_select[cols].description='Slope of the O/H gradient using the calibrator '+cal
        else:
            if (cols.find("e_OH")>-1):
                cal=cols.replace('e_OH_','')
                cal=cal.replace('_cen','')            
                tab_pe_SAMI_select[cols].description='Error in Oxygen abundance using the calibrator '+cal+' at the central region'
            else:
                cal=cols.replace('OH_','')
                cal=cal.replace('_cen','')            
                tab_pe_SAMI_select[cols].description='Oxygen abundance using the calibrator '+cal+' at the central region'
  
#print(tab_pe_SAMI_select['name'].description)
#print(tab_pe_SAMI_select['name'].units)
#print(tab_pe_SAMI_info['UNITS'])
