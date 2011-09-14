#!/usr/bin/env python
# encoding: utf-8
"""
interactive.py

Created by Jonathan Whitmore on 2010-06-14.
"""

from config_calibration import *

from enthought.traits.api \
    import HasTraits, Array, Range, Float, Enum, on_trait_change, Property
from enthought.traits.ui.api import View, Item
from enthought.chaco.chaco_plot_editor import ChacoPlotItem
from numpy import arange

# =============
# = Constants =
# =============
c_light = 299792458.


class Data(HasTraits):
    global initial_shift
    def initial_shift(fmultiple,fshift,fsigma):
        better_flx = flx * fmultiple
        better_kernel = normal_gaussian(elements,fsigma)
        better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                               slice_iof, mode='same'))
        better_y2 = si.splev(wav+fshift,better_tck)
        return sum((better_y2 - better_flx) ** 2 / \
                   (fmultiple * err) ** 2 )

    global class_shift
    def class_shift(fmultiple,fshift,fsigma):
       better_flx = starting_flx * fmultiple
       better_kernel = normal_gaussian(gauss_elements,fsigma)
       better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                              slice_iof, mode='same'))
       better_y2 = si.splev(wav+fshift,better_tck)
       return sum((better_y2 - better_flx) ** 2 / \
                  (fmultiple * err) ** 2 )
        

    def second_shift(fmultiple,fshift,fsigma,fskew):
        better_flx = starting_flx * fmultiple
        better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
        better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                               slice_iof, mode='same'))
        better_y2 = si.splev(wav+fshift,better_tck)
        return sum((better_y2 - better_flx) ** 2 / \
                   (fmultiple * err) ** 2 )

    # bin-shifting
    def shiftperbin(fshift,i,fmultiple,fsigma):
        better_flx = starting_flx[whit_index[int(i)]] * fmultiple
        better_kernel = normal_gaussian(gauss_elements,fsigma)
        better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                               slice_iof, mode='same'))
        better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
        return sum((better_y2 - better_flx) ** 2 / \
                    (fmultiple * err[whit_index[int(i)]]) ** 2)

    def skewshiftperbin(fshift,i,fmultiple,fsigma,fskew):
        better_flx = starting_flx[whit_index[int(i)]] * fmultiple
        better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
        better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                               slice_iof, mode='same'))
        better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
        return sum((better_y2 - better_flx) ** 2 / \
                    (fmultiple * err[whit_index[int(i)]]) ** 2)

    #
    def doshit(view):
        CalibrationSummaryDir = "Summary/calibration."

        # ====================
        # = Read in FTS Data =
        # ====================

        FTSReader = csv.reader(open(FTSFile), delimiter=' ')
        iow = []
        iof = []
        for row in FTSReader:
            iow.append(float(row[0]))
            iof.append(float(row[1]))

        iow = np.array(iow)
        iof = np.array(iof)
        today = datetime.date.today()
        CalibrationLogFile = "Logs/calibration." + str(today) + ".ascii"
        CalibrationLogFileHandle = open(CalibrationLogFile, 'a')    

        print "Beginning calibration analysis of data..."

        for core_exposure in exposures_to_analyze:
            CalibrationDataList = glob.glob(OutputContinuumDir + astro_object + "." + core_exposure + "*")
            print datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")

            for CalibrationDataFile in CalibrationDataList:
                #   Iterate over the Data files
                
                global wav
                global flx, err, pix
                wav = []
                flx = []
                err = []
                pix = []

                CalibrationDataReader = csv.reader(open(CalibrationDataFile), delimiter=' ')
                for row in CalibrationDataReader:
                    wav.append(float(row[0]))
                    flx.append(float(row[1]))
                    err.append(float(row[2]))
                    if len(row) > 3:
                        pix.append(float(row[3]))


                wav = np.array(wav)
                flx = np.array(flx)
                err = np.array(err)
                pix = np.array(pix)
                if len(wav) < 100:
                    continue

                # ================================================
                # = Logic Test to see if file overlaps FTS data =
                # ================================================
                if(wav[0] < iow[0]+10 or wav[-1] > iow[-1] - 100):
                    print CalibrationDataFile, "is outside of overlap"
                    continue

                # =================================
                # = Slice FTS to manageable size  =
                # =================================
                global slice_iow
                global slice_iof
                slice_iow = iow[np.where(iow > wav[0] - 10)]
                slice_iow = slice_iow[np.where(slice_iow < wav[-1] + 10)]
                slice_iof = iof[np.where(iow > wav[0] - 10)]
                slice_iof = slice_iof[np.where(slice_iow < wav[-1] + 10)]
                
                global starting_flx
                global elements
                elements = gauss_elements
                starting_flx = flx

                firstsplit = CalibrationDataFile.split("/")
                secondsplit = firstsplit[-1].split(".")
                naming = "." + secondsplit[1] + "." + secondsplit[2] + "." + secondsplit[3] 
                # print wav[0],wav[1]
                print test34
                m = mi.Minuit(view.class_shift(),fmultiple=0.82,\
                                            fshift=0.01,\
                                            fsigma=15.,\
                                            # elements=gauss_elements,\
                                            # fix_elements=True,\
                                            # tel_wav=wav,\
                                            # fix_tel_wav=True,\
                                            # tel_flx=flx,\
                                            # fix_tel_flx=True,\
                                            # slice_iow=slice_iow,\
                                            # fix_slice_iow=True,\
                                            # slice_iof=slice_iof,\
                                            # fix_slice_iof=True,\
                                            strategy=2)

                # m.printMode=1
                m.migrad()
                #m.minos()
                order_velocity = c_light / wav[len(wav)/2]
                print CalibrationDataFile, "\nShift: ", m.values["fshift"] * order_velocity, \
                        m.errors["fshift"] * order_velocity
                print "Sigma: ", m.values["fsigma"], m.errors["fsigma"]
        return 
    volume = Array
    pressure = Property(Array, depends_on=['temperature', 'attraction',
                                       'totVolume'])
    attraction = Range(low=-50.0,high=50.0,value=0.0)
    totVolume = Range(low=.01,high=100.0,value=0.01)
    temperature = Range(low=-50.0,high=50.0,value=50.0)
    r_constant= Float(8.314472)
    plot_type = Enum("line", "scatter")

    traits_view = View(ChacoPlotItem("volume", "pressure",
                               type_trait="plot_type",
                               resizable=True,
                               x_label="Volume",
                               y_label="Pressure",
                               x_bounds=(-10,120),
                               x_auto=False,
                               y_bounds=(-2000,4000),
                               y_auto=False,
                               color="blue",
                               bgcolor="white",
                               border_visible=True,
                               border_width=1,
                               title='Pressure vs. Volume',
                               padding_bg_color="lightgray"),
                       Item(name='attraction'),
                       Item(name='totVolume'),
                       Item(name='temperature'),
                       Item(name='r_constant', style='readonly'),
                       Item(name='plot_type'),
                       resizable = True,
                       buttons = ["OK"],
                       title='Interactive Calibration',
                       width=900, height=800)
    global test34
    test34 = 56
    

    def _volume_default(self):
        """ Default handler for volume Trait Array. """
        return arange(.1, 100)
        # return wav

    def _get_pressure(self):
        """Recalculate when one a trait the property depends on changes."""
        return ((self.r_constant*self.temperature)
              /(self.volume - self.totVolume)
             -(self.attraction/(self.volume*self.volume)))



#
        


if __name__ == '__main__':
    viewer = Data()
    viewer.configure_traits()
    viewer.doshit()