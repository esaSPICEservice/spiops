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
                           '$KERNELS/ck/ATNM_P170331235152_01351.BC'
                           '$KERNELS/ck/ATNM_VMC_170101_170712_V01.BC'

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
                           '$KERNELS/spk/ORMM_T19_170101000000_01318.BSP'
                           '$KERNELS/spk/ORMM_T19_170201000000_01326.BSP'
                           '$KERNELS/spk/ORMM_T19_170301000000_01336.BSP'
                           '$KERNELS/spk/ORMM_T19_170401000000_01344.BSP'
                           '$KERNELS/spk/ORMM_T19_170501000000_01351.BSP'
                           '$KERNELS/spk/ORMM_T19_170601000000_01351.BSP'

                         )

   \begintext


End of MK file.