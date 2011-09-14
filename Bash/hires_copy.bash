#!/bin/bash

# Directories to be set
RAW_DIR="Raw/"

for dirs in $RAW_DIR
do
    if [ -d $dirs ]; then
        echo "Directory $dirs exists..."
    else
        echo "Error: '$dirs' does not exist!!"
        echo "Put raw data into a Raw directory."
        echo "EXITING"
        exit 1
    fi
done

mkdir Info Logs 

mkdir Ascii_Output
mkdir Ascii_Output/Raw_Copy Ascii_Output/Raw_Masked
mkdir Ascii_Output/Continuum Ascii_Output/Calibration Ascii_Output/Resolution
mkdir Ascii_Output/Chip

mkdir PDF_Output
mkdir PDF_Output/Continuum PDF_Output/Calibration PDF_Output/Analysis

mkdir PS_Output
mkdir PS_Output/Continuum PS_Output/Calibration PS_Output/Analysis

mkdir QA_PDF
mkdir QA_PDF/Order-Residuals
mkdir QA_PDF/Order-Chi-Square
mkdir QA_PDF/Bins-Residuals
mkdir QA_PDF/Bins-Overlay
mkdir QA_PDF/Bins-Chi-Square

mkdir QA_PS
mkdir QA_PS/Order-Residuals
mkdir QA_PS/Order-Chi-Square
mkdir QA_PS/Bins-Residuals
mkdir QA_PS/Bins-Overlay
mkdir QA_PS/Bins-Chi-Square

mkdir Summary

mkdir Monte-Carlo

CLEAN_DIR="Ascii_Output/Raw_Copy/"

for dirs in $CLEAN_DIR
do
    if [ -d $dirs ]; then
        echo "Directory $dirs exists..."
    else
        echo "Error: '$dirs' does not exist!!"
        echo "Did you run Bash/setup.bash?"
        echo "EXITING"
        exit 2
    fi
done


echo "Copying Raw Data..."

for file in $(ls Raw/)
do
    # echo "$file" | awk -F _ '{print $2".k."$3}' | awk -F . '{print $1"."$2"."$3".ascii"}'
    # The first awk statement uses _ as the field separator. Make the $ correspond to exposure, chip, order.
    TEMP="$(echo "$file" | awk -F _ '{print $2".k."$3}' | awk -F . '{print $1"."$2"."$3".ascii"}')"
    # If Pixels are listed, print them fourth, otherwise, wavelength, flux, error, 
    awk '{if ( $NF == $3) print $1,$2,$3; if ( $NF == $4) print $2,$3,$4,$1}' "Raw/$file" >  $CLEAN_DIR$TEMP 
done

echo "Done"
