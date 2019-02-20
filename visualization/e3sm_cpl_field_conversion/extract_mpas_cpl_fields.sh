#!/bin/bash

usage()
{
cat<<EOF
usage: $0 options
This script:
1) extracts MPAS[LI,O,CICE]-generated/directed input/output coupler fields from coupler history files.
2) converts the dimension naming to MPAS standards.
3) appends all these fields, along with relevant mesh descriptor fields, to component-segregated output file(s) MPAS[LI,O,CICE]_CPL_fields.nc.
mpas*_CPL_fields.nc can then be read by, e.g., ParaView.
The user chooses which MPAS component fields (i.e. from MPASLI, MPASO, or MPASCICE) to extract by specifing one or more of the g,o,i flags (see below).  For example, if only -g is specified, the routine will *only* try to extract MPASLI-generated/directed input/output coupler fields.
-h show this message and exit
-c coupler file containing MPAS-based coupler fields (mandatory)
-g file containing MPAS-LI mesh descriptor information (optional)
-o file containing MPAS-O mesh descriptor information (optional)
-i file containing MPAS-CICE mesh descriptor information (optional)
Original author: Jeremy Fyke (fyke@lanl.gov)
Maintained: Matt Hoffman (mhoffman@lanl.gov)
EOF
}

MpasliOutFile=MPASLI_CPL_fields.nc
MpasOOutFile=MPASO_CPL_fields.nc
MpasCICEOutFile=MPASCICE_CPL_fields.nc

while getopts "hc:g:o:i" OPTION
do
   case $OPTION in
     h)
       usage
       exit 1
       ;;
     c)
       CPLFile=$OPTARG
       ;;
     g)
       MpasliGridFile=$OPTARG
       ;;
     o)
       MpasOGridFile=$OPTARG
       ;;
     i)
       MpasCICEGridFile=$OPTARG
       ;;
   esac
done

if [ -z $CPLFile ]; then
   echo Error: no coupler history file identified.
   exit 0
fi
 if [ -f $CPLFile ]; then
    echo ''
    echo Extracting fields from $CPLFile
 else
    echo Error: $CPLFile not found.
    exit 0
 fi

MPASMeshVars='angleEdge,areaCell,areaTriangle,cellsOnCell,cellsOnEdge,cellsOnVertex,dcEdge,dvEdge,edgesOnCell,edgesOnEdge,edgesOnVertex,indexToCellID,indexToEdgeID,indexToVertexID,kiteAreasOnVertex,latCell,latEdge,latVertex,lonCell,lonEdge,lonVertex,nEdgesOnCell,nEdgesOnEdge,verticesOnCell,verticesOnEdge,weightsOnEdge,xCell,xEdge,xVertex,yCell,yEdge,yVertex,zCell,zEdge,zVertex'

###MPAS-LI FIELDS###
if [ "$MpasliGridFile" ]; then
   if [ -f $MpasliGridFile ]; then
      echo ''
      echo Merging MPAS-LI-based coupler fields with MPAS-LI grid information from $MpasliGridFile.

      #Extract g2x coupler history fields to temporary file.
      ncks -O -v g2x+ $CPLFile g2xfile.nc
      ncrename -O -d g2x_nx,nCells g2xfile.nc g2xfile.nc
      ncrename -O -d time,Time g2xfile.nc g2xfile.nc
      ncwa -O -a g2x_ny g2xfile.nc g2xfile.nc

      #Extract x2g coupler history fields to temporary file.
      ncks -O -v x2g+ $CPLFile x2gfile.nc
      ncrename -O -d x2g_nx,nCells x2gfile.nc x2gfile.nc
      ncrename -O -d time,Time x2gfile.nc x2gfile.nc
      ncwa -O -a x2g_ny x2gfile.nc x2gfile.nc

      #Stick everything together in final output file.
      echo ''
      echo Outputting MPAS-LI-based coupler fields to $MpasliOutFile.

      #First, extract grid-descriptor fields from file
      ncks -O -v $MPASMeshVars $MpasliGridFile $MpasliOutFile

      #Then, tack on coupler fields
      ncks -A g2xfile.nc $MpasliOutFile
      ncks -A x2gfile.nc $MpasliOutFile
      ncap2 -s 'defdim("nVertLevels",1)' $MpasliOutFile       #remake vertical dimension (since not automatically obtained via strictly 2D coupler fields)

      #Clean up.
      rm g2xfile.nc x2gfile.nc
   else
      echo Error: $MpasliGridFile not found.
      exit 0
   fi
fi

###MPAS-O FIELDS###
if [ "$MpasOGridFile" ]; then
   if [ -f $MpasOGridFile ]; then
      echo ''
      echo Merging MPAS-O-based coupler fields with MPAS-O grid information from $MpasOGridFile.

      # Note: the . preceding dimensions in ncrename below suppresses an error if that dimension is not present
      # (which can be the case for certain compsets)

      #Extract o2x coupler history fields to temporary file.
      ncks -O -v o2x_+ $CPLFile o2xfile1.nc
      ncrename -O -d .o2x_nx,nCells o2xfile1.nc o2xfile1.nc
      ncrename -O -d time,Time o2xfile1.nc o2xfile1.nc
      ncwa -O -a o2x_ny o2xfile1.nc o2xfile1.nc

      #Extract o2xa coupler history fields to temporary file.
      ncks -O -v o2xa_+ $CPLFile o2xfile2.nc
      ncrename -O -d .o2xa_nx,nCells o2xfile2.nc o2xfile2.nc
      ncrename -O -d time,Time o2xfile2.nc o2xfile2.nc
      ncwa -O -a o2xa_ny o2xfile2.nc o2xfile2.nc

      #Extract x2oacc coupler history fields to temporary file.
      ncks -O -v x2oacc+ $CPLFile x2ofile.nc
      ncks -O -x -v x2oacc_ox_cnt x2ofile.nc x2ofile.nc #Prune counter variable
      ncrename -O -d .x2oacc_nx,nCells x2ofile.nc x2ofile.nc
      ncrename -O -d time,Time x2ofile.nc x2ofile.nc
      ncwa -O -a x2oacc_ny x2ofile.nc x2ofile.nc

      #Stick everything together in final output file.
      echo ''
      echo Outputting MPAS-O-based coupler fields to $MpasOOutFile.

      #First, extract grid-descriptor fields from file
      ncks -O -v $MPASMeshVars $MpasOGridFile $MpasOOutFile #extract grid-descriptor fields from file

      #Then, tack on coupler fields
      ncks -A o2xfile1.nc $MpasOOutFile
      ncks -A o2xfile2.nc $MpasOOutFile
      ncks -A x2ofile.nc $MpasOOutFile

      ncap2 -s 'defdim("nVertLevels",1)' $MpasOOutFile      #remake vertical dimension (since not automatically obtained via strictly 2D coupler fields)

      #Clean up.
      rm o2xfile1.nc o2xfile2.nc x2ofile.nc
   else
      echo Error: $MpasOGridFile not found.
      exit 0
   fi
fi

###MPAS-CICE FIELDS###
if [ "$MpasCICEGridFile" ]; then
   if [ -f $MpasCICEGridFile ]; then
      echo ''
      echo Merging MPAS-CICE-based coupler fields with MPAS-CICE grid information from $MpasCICEGridFile.
      echo ''
      echo WARNING: MPAS-CICE FIELD EXTRACTION NOT TESTED.
      echo ''

      #Extract i2x coupler history fields to temporary file.
      ncks -O -v i2x+ $CPLFile i2xfile.nc
      ncrename -O -d i2x_nx,nCells i2xfile.nc i2xfile.nc
      ncrename -O -d time,Time i2xfile.nc i2xfile.nc
      ncwa -O -a i2x_ny i2xfile.nc i2xfile.nc

      #Extract x2i coupler history fields to temporary file.
      ncks -O -v x2i+ $CPLFile x2ifile.nc
      ncrename -O -d x2i_nx,nCells x2ifile.nc x2ifile.nc
      ncrename -O -d time,Time x2ifile.nc x2ifile.nc
      ncwa -O -a x2i_ny x2ifile.nc x2ifile.nc

      #Stick everything together in final output file.
      echo ''
      echo Outputting MPAS-CICE-based coupler fields to $MpasCICEOutFile.

      #First, extract grid-descriptor fields from file
      ncks -O -v $MPASMeshVars $MpasCICEGridFile $MpasCICEOutFile #extract grid-descriptor fields from file

      #Then, tack on coupler fields
      ncks -A i2xfile.nc $MpasCICEOutFile
      ncks -A x2ifile.nc $MpasCICEOutFile

      ncap2 -s 'defdim("nVertLevels",1)' $MpasCICEOutFile         #remake vertical dimension (since not automatically obtained via strictly 2D coupler fields)

      #Clean up.
      rm i2xfile.nc x2ifile.nc
   else
      echo Error: $MpasCICEGridFile not found.
      exit 0
   fi
fi
