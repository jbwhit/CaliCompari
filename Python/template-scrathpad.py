



# Linear first

# x-shift 
d_x = 82.0

x1 = np.arange(100)
x2 = np.arange(100) + 1 * d_x
x3 = np.arange(100) + 2 * d_x


y1 = - (x1 -38) ** 2 + 250
y2 = - (x1 -38) ** 2 + 250
y3 = - (x1 -38) ** 2 + 250

xt = np.arange(100)
yt = - (x1 -38) ** 2 + 250



# delta y
d_y = 35.0

# m = mi.Minuit(initial_shift,fmultiple=0.82,fshift=0.01,\
                     # fsigma=15.,strategy=2)

  
def xshift(shift):
  """docstring for xshift"""
  return (np.average(x1) - np.average(xt - shift)) ** 2

mobject = mi.Minuit(xshift, shift = 45.0, strategy=2)

mobject.migrad()

def xshift2(shift):
  """docstring for xshift"""
  return (np.average(x2) - np.average(xt - shift)) ** 2

mobject2 = mi.Minuit(xshift2, shift = 45.0, strategy=2)
mobject2.migrad()
mobject2.values["shift"]

# m.printMode=1
# m.migrad()

XList = [x1,x2,x3]
YList = [y1,y2,y3]

for item in XList:
  print (np.average(item) - np.average(xt)) ** 2

List0064 = ['0064.k.058',
 '0064.k.059',
 '0064.k.060',
 '0064.k.061',
 '0064.k.062',
 '0064.k.063',
 '0064.k.064',
 '0064.k.065',
 '0064.k.066',
 '0064.k.067',
 '0064.k.068',
 '0064.k.069',
 '0064.k.070']

List0065 = ['0065.k.058',
  '0065.k.059',
  '0065.k.060',
  '0065.k.061',
  '0065.k.062',
  '0065.k.063',
  '0065.k.064',
  '0065.k.065',
  '0065.k.066',
  '0065.k.067',
  '0065.k.068',
  '0065.k.069',
  '0065.k.070']

List0068 = ['0068.k.058',
  '0068.k.059',
  '0068.k.060',
  '0068.k.061',
  '0068.k.062',
  '0068.k.063',
  '0068.k.064',
  '0068.k.065',
  '0068.k.066',
  '0068.k.067',
  '0068.k.068',
  '0068.k.069',
  '0068.k.070']

List0069 = ['0069.k.058',
  '0069.k.059',
  '0069.k.060',
  '0069.k.061',
  '0069.k.062',
  '0069.k.063',
  '0069.k.064',
  '0069.k.065',
  '0069.k.066',
  '0069.k.067',
  '0069.k.068',
  '0069.k.069',
  '0069.k.070']


deltax = []
multix = []
deltay = []
multiy = []

for item in List0064:
    deltax.append(np.average(avwav[item]) - np.average(xt))
    multix.append(avwav[item][-1] - avwav[item][0])
    deltay.append(np.average(avcal[item]))
    multiy.append(np.max(avcal[item] - np.min(avcal[item])))
    
  # mobject2 = mi.Minuit(xshift2, shift = 45.0, strategy = 2)
# 
# 
# avhist = {}
# avhistedge = {}
# avhistmid = {}
# avhistsize = {}
# bins = 10
# avhist[key1] = []
# key1hist = []
# key2hist = []
# 
# bin1edges = np.linspace(np.min(avwav[key1]),np.max(avwav[key1]),bins)
# bin2edges = np.linspace(np.min(avwav[key2]),np.max(avwav[key2]),bins)
# for i in range(bins-1):
#   bin1mid = (bin1edges[i + 1] + bin1edges[i])/2.0
#   bin2mid = (bin2edges[i + 1] + bin2edges[i])/2.0
#   bin1size = (bin1edges[i + 1] - bin1edges[i])/2.0
#   bin2size = (bin2edges[i + 1] - bin2edges[i])/2.0
#   key1hist.append(np.average(avcal[key1][np.where(np.abs(avwav[key1] - bin1mid) - bin1size < 0)]))
#   key2hist.append(np.average(avcal[key2][np.where(np.abs(avwav[key2] - bin2mid) - bin2size < 0)]))

# scaling histogram proportionally 
avhist = {}
avhistedge = {}
avhistmid = {}
avhistsize = {}
bins = 10
for key in List0064:
  avhist[key] = []
  avhistedge[key] = np.linspace(np.min(avwav[key]),np.max(avwav[key]),bins)
  for i in range(bins-1):
    avhistmid[key] = (avhistedge[key][i + 1] + avhistedge[key][i]) / 2.0
    avhistsize[key] = (avhistedge[key][i + 1] - avhistedge[key][i]) / 2.0
    avhist[key].append(np.average(avcal[key][np.where(np.abs(avwav[key] - avhistmid[key] ) - avhistsize[key] < 0)]))


# Keeping wavelength stuff similar
# Error when running over the second item in the list
avhist = {}
avhistedge = {}
avhistmid = {}
avhistsize = {}
bins = 15
avhistedge['0064.k.058'] = np.linspace(np.min(avwav['0064.k.058'][0] - np.average(avwav['0064.k.058'])),np.max(avwav['0064.k.058'][-1] - np.average(avwav['0064.k.058'])),bins)
for key in List0064:
  avhist[key] = []
  for i in range(bins-1):
    avhistmid[key] = (avhistedge['0064.k.058'][i + 1] + avhistedge['0064.k.058'][i]) / 2.0
    avhistsize[key] = (avhistedge['0064.k.058'][i + 1] - avhistedge['0064.k.058'][i]) / 2.0
    avhist[key].append(np.average(avcal[key][np.where(np.abs(avwav[key] -np.average(avwav[key]) - avhistmid[key]) - avhistsize[key] < 0)]))

tempx = []
tempy = []
for key in List0064:
  tempx.append(np.average(avwav[key]))
  tempy.append(np.max(avcal[key]))

regresstuple = linregress(tempx,tempy)
for key in List0064:
  pl.plot(avhist[key] - np.average(avwav[key]) * regresstuple[0] - regresstuple[1])

mx2 = []
my2 = []
for key in List0064:
  mx2.append(np.average(avwav[key]))
  my2.append(np.max(avcal[key]) - np.min(avcal[key]))

# I'm dividing the calibration by the change in height after bringing everything level
mtuple = linregress(mx2,my2)
for key in List0064:
  pl.plot((avcal[key] - np.average(avwav[key]) * regresstuple[0] - regresstuple[1])/(np.average(avwav[key] * mtuple[0] + mtuple[1])))

# The number of points in the lowest order happens to be the number of points in all the orders.
template = []
for i in range(len(avcal['0064.k.058'])):
  hold = []
  for key in List0064:
    hold.append((avcal[key][i] - np.average(avwav[key]) * regresstuple[0] - regresstuple[1])/(np.average(avwav[key] * mtuple[0] + mtuple[1])))
  template.append(np.average(hold))

plate = np.array(template)
for key in List0064:
  pl.plot(avwav[key], plate * (np.average(avwav[key]) * mtuple[0] + mtuple[1] )  + np.average(avwav[key]) * regresstuple[0] + regresstuple[1])
  pl.plot(avwav[key],avcal[key])

# Besides the averages, there are 4 numbers, the slope and y-intercept for 2 linear fits: 
# 1. The line through the maximum vshift as a function of average wavelength in an echelle.
# 2. The line through the difference between the maximum and minimum calibration values as a function of average wavelength in an echelle.

tempx = []
tempy = []
for key in List0064:
  tempx.append(np.average(avwav[key]))
  tempy.append(np.max(avcal[key]))

# Slope and y-intercept of max calibration value and average wavelength
var1, var2, trash, trash, trash = linregress(tempx,tempy)

mx2 = []
my2 = []
for key in List0064:
  mx2.append(np.average(avwav[key]))
  my2.append(np.max(avcal[key]) - np.min(avcal[key]))

# Slope and y-intercept of overall mulitplication factor 
var3, var4, trash, trash, trash = linregress(mx2,my2)

# The number of points in the lowest order happens to be the number of points in all the orders.

def sumsqdiff(maxslope, maxint, multslope, multint, FitExposureArray, TestExposureArray):
  """docstring for sumsqdiff"""
  sqdiff = []
  plate = template_fitting(maxslope, maxint, multslope, multint, FitExposureArray)
  for key in TestExposureArray:
    sqdiff.append(np.sum((avcal[key] - (plate * (np.average(avwav[key]) * multslope + multint )  + np.average(avwav[key]) * maxslope + multint))**2))
  return np.sum(sqdiff)

def error_function(v, FitExposureArray, TestExposureArray):
  """docstring for error_function"""
  sqdiff = []
  plate = template_fitting(v[0], v[1], v[2], v[3], FitExposureArray)
  for key in TestExposureArray:
    sqdiff.append(np.sum((avcal[key] - (plate * (np.average(avwav[key]) * v[2] + v[3] ) + np.average(avwav[key]) * v[0] + v[1]))**2))
  return np.sum(sqdiff)

def error_function_array(v, FitExposureArray, TestExposureArray):
  """docstring for error_function"""
  sqdiff = []
  plate = template_fitting(v, FitExposureArray)
  for key in TestExposureArray:
    sqdiff.append(avcal[key] - (plate * (np.average(avwav[key]) * v[2] + v[3] ) - (np.average(avwav[key]) * v[0] + v[1])))
  return np.concatenate(sqdiff)

def template_fitting(v, FitExposureArray):
  """Fitting template for Keck runs"""
  template = []
  for i in range(len(avcal['0064.k.058'])):
    hold = []
    for key in FitExposureArray:
      hold.append((avcal[key][i] - np.average(avwav[key]) * v[0] - v[1])/(np.average(avwav[key] * v[2] + v[3])))
    template.append(np.average(hold))  
  return np.array(template)

def plotfit(v, FitExposureArray, TestExposureArray):
  """docstring for plotfit"""
  plate = template_fitting(v, FitExposureArray)
  for key in TestExposureArray:
    pl.plot(avwav[key], (plate * (np.average(avwav[key]) * v[2] + v[3] ) - (np.average(avwav[key]) * v[0] + v[1])))
    pl.plot(avwav[key],avcal[key])
  pass

def plotresiduals(v, FitExposureArray, TestExposureArray):
  """docstring for plotresiduals"""
  plate = template_fitting(v, FitExposureArray)
  for key in TestExposureArray:
    pl.plot(avwav[key], avcal[key] - (plate * (np.average(avwav[key]) * v[2] + v[3] ) - (np.average(avwav[key]) * v[0] + v[1])))
  pass

def linearresiduals(FitExposureArray):
  """docstring for linearresiduals"""
  temporaryx = []
  temporaryy = []
  for key in FitExposureArray:
    temporaryx.append(avwav[key])
    temporaryy.append(avcal[key])
    
  temporaryx = np.concatenate(temporaryx)
  temporaryy = np.concatenate(temporaryy)
  print temporaryx
  slope = linregress(temporaryx,temporaryy)[0]
  intercept = linregress(temporaryx,temporaryy)[1]
  pl.plot(temporaryx, temporaryx * slope + intercept)
  pass

def explainv(v):
  """docstring for explainv"""
  print 
  pass
  
  
# Example of how to fit least squares to function
# v0 = 
from scipy.optimize import leastsq

v0065, success = leastsq(error_function_array, v0, args=(List0065, List0065), maxfev=10000)

def quadratic_error_function_array(v, FitExposureArray, TestExposureArray):
  """docstring for error_function"""
  sqdiff = []
  plate = template_fitting(v, FitExposureArray)
  for key in TestExposureArray:
    sqdiff.append(avcal[key] - (plate * (np.average(avwav[key]) * v[2] + v[3] ) - (np.average(avwav[key]) * v[0] + v[1])))
  return np.concatenate(sqdiff)
    
template = []
for i in range(len(avcal['0064.k.058'])):
  hold = []
  for key in List0064:
    hold.append((avcal[key][i] - np.average(avwav[key]) * regresstuple[0] - regresstuple[1])/(np.average(avwav[key] * mtuple[0] + mtuple[1])))
  template.append(np.average(hold))

plate = np.array(template)
sqdiff = []
for key in List0064:
  sqdiff.append(np.sum((avcal[key] - (plate * (np.average(avwav[key]) * mtuple[0] + mtuple[1] )  + np.average(avwav[key]) * regresstuple[0] + regresstuple[1]))**2))

sumsqdiff = np.sum(sqdiff)

