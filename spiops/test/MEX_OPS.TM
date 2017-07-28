KPL/MK

Meta-kernel for Mars Express dataset v0.2.0 -- Operational 20170613_001
==========================================================================

   This meta-kernel lists the Mars Express Operational SPICE kernels
   providing information for the full mission based on predicted, test
   and/or measured data. This set of kernels and the order in which
   they are listed are picked to provide the best available data and
   the most complete coverage based on the information about the
   kernels available at the time this meta-kernel was made for the
   Mars Express Operational scenario.

     
Usage of the Meta-kernel                               
-------------------------------------------------------------------------    

   The SPK files present in this meta-kernel have the following coverage:

         Body: MARS EXPRESS (-41)
            Start of Interval (TDB)             End of Interval (TDB)
            -----------------------------       -----------------------------
            2003-06-02T19:12:20.596             2017-06-30T20:47:50.179


   The CK files present in this meta-kernel have the following coverage:

         Body: MEX_SC_REF (-41001)
            Start of Interval (TDB)             End of Interval (TDB)
            -----------------------------       -----------------------------
            2006-04-01T00:00:00.000             2017-06-27T06:09:29.157


   Reading the comments of the binary SPK and C-Kernels is very helpful to 
   understand the implemented S/C position and orientation. You can use the
   following NAIF command line utility to extract that information from the 
   kernels:

      > commnt -r MEX_*

   The kernels listed below can be obtained from the ESA SPICE FTP server:
      
      ftp://spiftp.esac.esa.int/data/SPICE/MARS-EXPRESS/kernels/


Implementation Notes
-------------------------------------------------------------------------    

   It is recommended that users make a local copy of this file and
   modify the value of the PATH_VALUES keyword to point to the actual
   location of the Mars Express SPICE data set's ``data'' directory on
   their system. Replacing ``/'' with ``\'' and converting line 
   terminators to the format native to the user's system may also be 
   required if this meta-kernel is to be used on a non-UNIX workstation.


-------------------   

   This file was created on June 13, 2017 by Marc Costa Sitja ESAC/ESA.


   \begindata

     PATH_VALUES       = ( '/Users/mcosta/Dropbox/SPICE/SPICE_MEX/ftp/data/SPICE/MARS-EXPRESS/kernels' )

     PATH_SYMBOLS      = ( 'KERNELS' )

     KERNELS_TO_LOAD   = (
                           '$KERNELS/ck/ATNM_P060401000000_01346.BC'
                           '$KERNELS/ck/ATNM_P170331235152_01351.BC'
                           '$KERNELS/ck/ATNM_RECONSTITUTED_00004.BC'
                           '$KERNELS/ck/MEX_ASPERA_SAF_030603_68.BC'
                           '$KERNELS/ck/MEX_SA_2003_V0001.BC'
                           '$KERNELS/ck/MEX_SA_2004_V0001.BC'
                           '$KERNELS/ck/MEX_SA_2005_V0001.BC'
                           '$KERNELS/ck/MEX_SA_2006_V0001.BC'
                           '$KERNELS/ck/MEX_SA_2007_V0001.BC'
                           '$KERNELS/ck/MEX_SA_2008_V0033.BC'
                           '$KERNELS/ck/MEX_SA_2009_V0054.BC'
                           '$KERNELS/ck/MEX_SA_2010_V0048.BC'
                           '$KERNELS/ck/MEX_SA_2011_V0048.BC'
                           '$KERNELS/ck/MEX_SA_2012_V0046.BC'
                           '$KERNELS/ck/MEX_SA_2013_V0048.BC'
                           '$KERNELS/ck/MEX_SA_2014_V0044.BC'
                           '$KERNELS/ck/MEX_SA_2015_V0002.BC'
                           '$KERNELS/ck/MEX_SA_2016_V0022.BC'

                           '$KERNELS/fk/EARTHFIXEDIAU.TF'
                           '$KERNELS/fk/EARTH_TOPO_050714.TF'
                           '$KERNELS/fk/ESTRACK_V01.TF'
                           '$KERNELS/fk/MEX_V12.TF'
                           '$KERNELS/fk/NEW_NORCIA_TOPO.TF'
                           '$KERNELS/fk/RSSD0002.TF'

                           '$KERNELS/ik/MEX_ASPERA_V06.TI'
                           '$KERNELS/ik/MEX_HRSC_V08.TI'
                           '$KERNELS/ik/MEX_OMEGA_V01.TI'
                           '$KERNELS/ik/MEX_PFS_V03.TI'
                           '$KERNELS/ik/MEX_VMC_V01.TI'

                           '$KERNELS/lsk/NAIF0012.TLS'

                           '$KERNELS/pck/DE403-MASSES.TPC'
                           '$KERNELS/pck/PCK00010.TPC'

                           '$KERNELS/sclk/MEX_170601_STEP.TSC'

                           '$KERNELS/spk/DE405.BSP'
                           '$KERNELS/spk/EARTHSTNS_ITRF93_050714.BSP'
                           '$KERNELS/spk/ESTRACK_V01.BSP'
                           '$KERNELS/spk/MAR097.BSP'
                           '$KERNELS/spk/MEX_ASPERA_STRUCT_V12.BSP'
                           '$KERNELS/spk/NEW_NORCIA.BSP'
                           '$KERNELS/spk/ORHM_______________00038.BSP'
                           '$KERNELS/spk/ORMM_T19_031222180906_00052.BSP'
                           '$KERNELS/spk/ORMM_T19_040201000000_00060.BSP'
                           '$KERNELS/spk/ORMM_T19_040301000000_00068.BSP'
                           '$KERNELS/spk/ORMM_T19_040401000000_00072.BSP'
                           '$KERNELS/spk/ORMM_T19_040501000000_00076.BSP'
                           '$KERNELS/spk/ORMM_T19_040601000000_00080.BSP'
                           '$KERNELS/spk/ORMM_T19_040701000000_00086.BSP'
                           '$KERNELS/spk/ORMM_T19_040801000000_00088.BSP'
                           '$KERNELS/spk/ORMM_T19_040901000000_00091.BSP'
                           '$KERNELS/spk/ORMM_T19_041001000000_00096.BSP'
                           '$KERNELS/spk/ORMM_T19_041101000000_00100.BSP'
                           '$KERNELS/spk/ORMM_T19_041201000000_00105.BSP'
                           '$KERNELS/spk/ORMM_T19_050101000000_00109.BSP'
                           '$KERNELS/spk/ORMM_T19_050201000000_00113.BSP'
                           '$KERNELS/spk/ORMM_T19_050301000000_00117.BSP'
                           '$KERNELS/spk/ORMM_T19_050401000000_00122.BSP'
                           '$KERNELS/spk/ORMM_T19_050501000000_00129.BSP'
                           '$KERNELS/spk/ORMM_T19_050601000000_00136.BSP'
                           '$KERNELS/spk/ORMM_T19_050701000000_00140.BSP'
                           '$KERNELS/spk/ORMM_T19_050801000000_00153.BSP'
                           '$KERNELS/spk/ORMM_T19_050901000000_00165.BSP'
                           '$KERNELS/spk/ORMM_T19_051001000000_00181.BSP'
                           '$KERNELS/spk/ORMM_T19_051101000000_00194.BSP'
                           '$KERNELS/spk/ORMM_T19_051201000000_00203.BSP'
                           '$KERNELS/spk/ORMM_T19_060101000000_00217.BSP'
                           '$KERNELS/spk/ORMM_T19_060201000000_00229.BSP'
                           '$KERNELS/spk/ORMM_T19_060301000000_00241.BSP'
                           '$KERNELS/spk/ORMM_T19_060401000000_00253.BSP'
                           '$KERNELS/spk/ORMM_T19_060501000000_00263.BSP'
                           '$KERNELS/spk/ORMM_T19_060601000000_00274.BSP'
                           '$KERNELS/spk/ORMM_T19_060701000000_00290.BSP'
                           '$KERNELS/spk/ORMM_T19_060801000000_00294.BSP'
                           '$KERNELS/spk/ORMM_T19_060901000000_00307.BSP'
                           '$KERNELS/spk/ORMM_T19_061001000000_00310.BSP'
                           '$KERNELS/spk/ORMM_T19_061101000000_00324.BSP'
                           '$KERNELS/spk/ORMM_T19_061201000000_00336.BSP'
                           '$KERNELS/spk/ORMM_T19_070101000000_00349.BSP'
                           '$KERNELS/spk/ORMM_T19_070201000000_00361.BSP'
                           '$KERNELS/spk/ORMM_T19_070301000000_00371.BSP'
                           '$KERNELS/spk/ORMM_T19_070401000000_00387.BSP'
                           '$KERNELS/spk/ORMM_T19_070501000000_00395.BSP'
                           '$KERNELS/spk/ORMM_T19_070601000000_00403.BSP'
                           '$KERNELS/spk/ORMM_T19_070701000000_00414.BSP'
                           '$KERNELS/spk/ORMM_T19_070801000000_00423.BSP'
                           '$KERNELS/spk/ORMM_T19_070901000000_00431.BSP'
                           '$KERNELS/spk/ORMM_T19_071001000000_00441.BSP'
                           '$KERNELS/spk/ORMM_T19_071101000000_00449.BSP'
                           '$KERNELS/spk/ORMM_T19_071201000000_00457.BSP'
                           '$KERNELS/spk/ORMM_T19_080101000000_00465.BSP'
                           '$KERNELS/spk/ORMM_T19_080201000000_00474.BSP'
                           '$KERNELS/spk/ORMM_T19_080301000000_00484.BSP'
                           '$KERNELS/spk/ORMM_T19_080401000000_00490.BSP'
                           '$KERNELS/spk/ORMM_T19_080501000000_00498.BSP'
                           '$KERNELS/spk/ORMM_T19_080601000000_00508.BSP'
                           '$KERNELS/spk/ORMM_T19_080701000000_00516.BSP'
                           '$KERNELS/spk/ORMM_T19_080801000000_00525.BSP'
                           '$KERNELS/spk/ORMM_T19_080901000000_00536.BSP'
                           '$KERNELS/spk/ORMM_T19_081001000000_00544.BSP'
                           '$KERNELS/spk/ORMM_T19_081101000000_00549.BSP'
                           '$KERNELS/spk/ORMM_T19_081201000000_00553.BSP'
                           '$KERNELS/spk/ORMM_T19_090101000000_00562.BSP'
                           '$KERNELS/spk/ORMM_T19_090201000000_00570.BSP'
                           '$KERNELS/spk/ORMM_T19_090301000000_00580.BSP'
                           '$KERNELS/spk/ORMM_T19_090401000000_00586.BSP'
                           '$KERNELS/spk/ORMM_T19_090501000000_00594.BSP'
                           '$KERNELS/spk/ORMM_T19_090601000000_00602.BSP'
                           '$KERNELS/spk/ORMM_T19_090701000000_00611.BSP'
                           '$KERNELS/spk/ORMM_T19_090801000000_00619.BSP'
                           '$KERNELS/spk/ORMM_T19_090901000000_00628.BSP'
                           '$KERNELS/spk/ORMM_T19_091001000000_00637.BSP'
                           '$KERNELS/spk/ORMM_T19_091101000000_00648.BSP'
                           '$KERNELS/spk/ORMM_T19_091201000000_00654.BSP'
                           '$KERNELS/spk/ORMM_T19_100101000000_00662.BSP'
                           '$KERNELS/spk/ORMM_T19_100201000000_00674.BSP'
                           '$KERNELS/spk/ORMM_T19_100301000000_00682.BSP'
                           '$KERNELS/spk/ORMM_T19_100401000000_00688.BSP'
                           '$KERNELS/spk/ORMM_T19_100501000000_00696.BSP'
                           '$KERNELS/spk/ORMM_T19_100601000000_00704.BSP'
                           '$KERNELS/spk/ORMM_T19_100701000000_00712.BSP'
                           '$KERNELS/spk/ORMM_T19_100801000000_00723.BSP'
                           '$KERNELS/spk/ORMM_T19_100901000000_00731.BSP'
                           '$KERNELS/spk/ORMM_T19_101001000000_00739.BSP'
                           '$KERNELS/spk/ORMM_T19_101101000000_00749.BSP'
                           '$KERNELS/spk/ORMM_T19_101201000000_00755.BSP'
                           '$KERNELS/spk/ORMM_T19_110101000000_00758.BSP'
                           '$KERNELS/spk/ORMM_T19_110201000000_00766.BSP'
                           '$KERNELS/spk/ORMM_T19_110301000000_00774.BSP'
                           '$KERNELS/spk/ORMM_T19_110401000000_00780.BSP'
                           '$KERNELS/spk/ORMM_T19_110501000000_00788.BSP'
                           '$KERNELS/spk/ORMM_T19_110601000000_00794.BSP'
                           '$KERNELS/spk/ORMM_T19_110701000000_00802.BSP'
                           '$KERNELS/spk/ORMM_T19_110801000000_00815.BSP'
                           '$KERNELS/spk/ORMM_T19_110901000000_00824.BSP'
                           '$KERNELS/spk/ORMM_T19_111001000000_00834.BSP'
                           '$KERNELS/spk/ORMM_T19_111101000000_00842.BSP'
                           '$KERNELS/spk/ORMM_T19_111201000000_00848.BSP'
                           '$KERNELS/spk/ORMM_T19_120101000000_00859.BSP'
                           '$KERNELS/spk/ORMM_T19_120201000000_00867.BSP'
                           '$KERNELS/spk/ORMM_T19_120301000000_00873.BSP'
                           '$KERNELS/spk/ORMM_T19_120401000000_00879.BSP'
                           '$KERNELS/spk/ORMM_T19_120501000000_00887.BSP'
                           '$KERNELS/spk/ORMM_T19_120601000000_00896.BSP'
                           '$KERNELS/spk/ORMM_T19_120701000000_00906.BSP'
                           '$KERNELS/spk/ORMM_T19_120801000000_00914.BSP'
                           '$KERNELS/spk/ORMM_T19_120901000000_00920.BSP'
                           '$KERNELS/spk/ORMM_T19_121001000000_00930.BSP'
                           '$KERNELS/spk/ORMM_T19_121101000000_00936.BSP'
                           '$KERNELS/spk/ORMM_T19_121201000000_00944.BSP'
                           '$KERNELS/spk/ORMM_T19_130101000000_00952.BSP'
                           '$KERNELS/spk/ORMM_T19_130201000000_00961.BSP'
                           '$KERNELS/spk/ORMM_T19_130301000000_00964.BSP'
                           '$KERNELS/spk/ORMM_T19_130401000000_00970.BSP'
                           '$KERNELS/spk/ORMM_T19_130501000000_00977.BSP'
                           '$KERNELS/spk/ORMM_T19_130601000000_00985.BSP'
                           '$KERNELS/spk/ORMM_T19_130701000000_00995.BSP'
                           '$KERNELS/spk/ORMM_T19_130801000000_01003.BSP'
                           '$KERNELS/spk/ORMM_T19_130901000000_01011.BSP'
                           '$KERNELS/spk/ORMM_T19_131001000000_01019.BSP'
                           '$KERNELS/spk/ORMM_T19_131101000000_01025.BSP'
                           '$KERNELS/spk/ORMM_T19_131201000000_01033.BSP'
                           '$KERNELS/spk/ORMM_T19_140101000000_01041.BSP'
                           '$KERNELS/spk/ORMM_T19_140201000000_01050.BSP'
                           '$KERNELS/spk/ORMM_T19_140301000000_01058.BSP'
                           '$KERNELS/spk/ORMM_T19_140401000000_01064.BSP'
                           '$KERNELS/spk/ORMM_T19_140501000000_01070.BSP'
                           '$KERNELS/spk/ORMM_T19_140601000000_01076.BSP'
                           '$KERNELS/spk/ORMM_T19_140701000000_01084.BSP'
                           '$KERNELS/spk/ORMM_T19_140801000000_01092.BSP'
                           '$KERNELS/spk/ORMM_T19_140901000000_01102.BSP'
                           '$KERNELS/spk/ORMM_T19_141001000000_01110.BSP'
                           '$KERNELS/spk/ORMM_T19_141101000000_01118.BSP'
                           '$KERNELS/spk/ORMM_T19_141201000000_01124.BSP'
                           '$KERNELS/spk/ORMM_T19_150101000000_01132.BSP'
                           '$KERNELS/spk/ORMM_T19_150201000000_01140.BSP'
                           '$KERNELS/spk/ORMM_T19_150301000000_01148.BSP'
                           '$KERNELS/spk/ORMM_T19_150401000000_01152.BSP'
                           '$KERNELS/spk/ORMM_T19_150501000000_01154.BSP'
                           '$KERNELS/spk/ORMM_T19_150601000000_01162.BSP'
                           '$KERNELS/spk/ORMM_T19_150701000000_01170.BSP'
                           '$KERNELS/spk/ORMM_T19_150801000000_01181.BSP'
                           '$KERNELS/spk/ORMM_T19_150901000000_01191.BSP'
                           '$KERNELS/spk/ORMM_T19_151001000000_01201.BSP'
                           '$KERNELS/spk/ORMM_T19_151101000000_01211.BSP'
                           '$KERNELS/spk/ORMM_T19_151201000000_01215.BSP'
                           '$KERNELS/spk/ORMM_T19_160101000000_01222.BSP'
                           '$KERNELS/spk/ORMM_T19_160201000000_01231.BSP'
                           '$KERNELS/spk/ORMM_T19_160301000000_01239.BSP'
                           '$KERNELS/spk/ORMM_T19_160401000000_01245.BSP'
                           '$KERNELS/spk/ORMM_T19_160501000000_01253.BSP'
                           '$KERNELS/spk/ORMM_T19_160601000000_01259.BSP'
                           '$KERNELS/spk/ORMM_T19_160701000000_01267.BSP'
                           '$KERNELS/spk/ORMM_T19_160801000000_01278.BSP'
                           '$KERNELS/spk/ORMM_T19_160901000000_01285.BSP'
                           '$KERNELS/spk/ORMM_T19_161001000000_01294.BSP'
                           '$KERNELS/spk/ORMM_T19_161101000000_01303.BSP'
                           '$KERNELS/spk/ORMM_T19_161201000000_01309.BSP'
                           '$KERNELS/spk/ORMM_T19_170101000000_01318.BSP'
                           '$KERNELS/spk/ORMM_T19_170201000000_01326.BSP'
                           '$KERNELS/spk/ORMM_T19_170301000000_01336.BSP'
                           '$KERNELS/spk/ORMM_T19_170401000000_01344.BSP'
                           '$KERNELS/spk/ORMM_T19_170501000000_01351.BSP'
                           '$KERNELS/spk/ORMM_T19_170601000000_01351.BSP'

                         )

   \begintext


End of MK file.