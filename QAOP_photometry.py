from photutils.aperture import aperture_photometry #main photometry package/method - inherently linked to the apertures we import later
from photutils.aperture import ApertureStats #For getting the median and area that we use to calculate the background
from astropy.table import Table #For handling tables
from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus
from astropy.coordinates import SkyCoord #for defining aperture positions
import astropy.units as u #Need 'deg' and 'arcsec' for skycoord and aperture
from astropy.io import fits #for loading fits files
from astropy.wcs import WCS #for getting WCS data from header
import os #for file handling and saving!

#------------------
#These paths should be changed if needed
ApertureFilePath = "apertures.csv"
PhotometryFilePath = "photometry_master.ecsv"
#------------------


def loadAperturesFromFile(ApertureFilePath="apertures.csv"):
    '''Creates a list of both apertures and annuli, as well as a list of names that they were assigned. It reads in from a file, which should be a table written by the STARID & RADIALPROF scripts in the QAOP package. Namely, this should be a table with the columns 
    Name | RA | DEC | r | r_in | r_out
    Where each line is a star location, with the name it was assigned, it's location in RA and DEC (ICRS), and as well, the radius of the stars aperture itself, as well as the inner and outer radius for the annulus that is used to determine the local background for that star.
    '''
    apertureTab = Table.read(ApertureFilePath)

    names,apertures,annuli = [],[],[]

    for row in apertureTab:
        #print(row)
        names.append(row['Name'])
        location = SkyCoord(ra=row['RA']*u.deg,dec=row['DEC']*u.deg)
        #print(location)
        aperture = SkyCircularAperture(location,r=row['r']*u.arcsec)
        #print(aperture)
        apertures.append(aperture)
        annulus = SkyCircularAnnulus(location,r_in=row['r_in']*u.arcsec,r_out=row['r_out']*u.arcsec)
        #print(annulus)
        annuli.append(annulus)
    return names, apertures, annuli

#The following 5 functions each perform a base component of the total photometry process, and are then wrapped into a single callable function in the 6th function.

def getRawSum(image,aperture,wcs):
    photResults = aperture_photometry(image,aperture,wcs=wcs)
    aperture_raw_sum = photResults["aperture_sum"].data[0] 
    #its ready for multiple apertures, so they're returned in an array
    #We've only got one item though
    return aperture_raw_sum

def getArea(image,aperture,wcs):
    #aperture_area = ApertureStats(image,aperture,wcs=wcs)["area"].data
    aperture_area = ApertureStats(image,aperture,wcs=wcs).sum_aper_area.value 
    #    we don't care that it's pix^2
    #aperture_area = aperture.to_pixel(wcs).area
    #print(aperture_area.value)
    return aperture_area

def getMedian(image,annulus,wcs):
    annulus_stats = ApertureStats(image,annulus,wcs=wcs)
    #aperture_median_background = annulus_stats["median"].data
    aperture_median_background = annulus_stats.median
    return aperture_median_background

def calcBackground(area,median):
    backgroud_to_subtract = area*median
    return backgroud_to_subtract

def subBackground(raw_sum, back_to_sub):
    adjusted_sum = raw_sum - back_to_sub
    return adjusted_sum

def photValWrapper(image,aperture,annulus,wcs):
    '''This function packages together all the base component info gathering into a single function, and will return a dictionary with the values for that aperture. This result is then added to a new dictionary, which maps the aperture names to the result dictionaries for each of the apertures. This is then saved to a master 3D table as a single row for the timestamp of the image.'''
    aperture_raw_sum = getRawSum(image,aperture,wcs)
    aperture_area = getArea(image,aperture,wcs)
    annulus_median = getMedian(image,annulus,wcs)
    aperture_background = calcBackground(aperture_area,annulus_median)
    aperture_sum = subBackground(aperture_raw_sum,aperture_background)
    resultDict = {}
    resultDict["aperture_raw_sum"] = aperture_raw_sum
    resultDict["aperture_area"] = aperture_area
    resultDict["annulus_median"] = annulus_median
    resultDict["background_to_subtract"] = aperture_background
    resultDict["aperture_sum"] = aperture_sum
    return resultDict

def doForApertures(image,names,apertures,annuli,wcs):
    '''A function that calls photValWrapper() for each of the apertures it is given, and saves the result dictionary returned from that call to a new dictionary where the key is the name assigned to the aperture it passed.
    
    The parameter "apertures" should be the result of the aperture preparation, namely, it should be a list of SkyCircularAperture in the same order as names and annuli, such as that returned by loadAperturesFromFile().'''
    image_results = {}
    for i,name in enumerate(names): #I figure I'll need the i to acces the things about the aperture I've currently got referenced by name
        aperture_name = name
        aperture_results = photValWrapper(image,apertures[i],annuli[i],wcs=wcs)
        image_results[aperture_name] = aperture_results #save the result dict to a dict with what it was the result for
    return image_results

#The next section of the program deals with File IO and thusly iterating through a selection of images and doing the photometry on each of them.

def loadImageAndWCS(filepath):
    with fits.open(filepath) as hdul:
        image = hdul[0].data
        wcs = WCS(hdul[0].header)
        img_time = hdul[0].header['date-obs']
        return image,wcs,img_time
    
def doForFile(filepath,names,apertures,annuli):
    '''This function is designed to be used and called by a wrapper iterating though a subset of files. '''
    image,wcs,img_time = loadImageAndWCS(filepath)
    #I need to do something better with the output from this; like saving it to a table
    return img_time, doForApertures(image,names,apertures,annuli,wcs)

class photInstance:#(apertureFilePath='apertures.csv',resultDir='photometry'):
    #names = [] #I believe these C# style declarations are automatically handled by __init__.
    #apertures = []
    #annuli = []
    master_history = []
    
    def __init__(self,apertureFilePath='apertures.csv',resultDir='photometry'):
        self.names,self.apertures,self.annuli = loadAperturesFromFile(apertureFilePath)
        #now, check there's a results dir in the root, and make one if not
        if not os.listdir().count(resultDir): os.mkdir(resultDir)
        #now that we know the directoru exists, we can safely check how many master_table are in it :)
        master_count = os.listdir(resultDir).count('master_table.ecsv')
        if master_count:
            self.master_tab = Table.read('photometry/master_table.ecsv')
        elif not master_count:
            self.master_tab = self.createMasterTable()
    
    def createMasterTable(self):
        MasterResultTab = Table(names=np.hstack(('Time',self.names)),dtype=np.hstack((str,np.full(len(self.names),dict))))
        MasterResultTab.write(self.resultDir+'/master_table.ecsv',format='ascii.ecsv')
        return MasterResultTab
    
    def runForFile(self,filepath):
        image,wcs,img_time = loadImageAndWCS(filepath)
        resultDict = doForApertures(image,self.names,self.apertures,self.annuli,wcs)
        resultDict['Time'] = img_time
        self.addRowToMaster(resultDict,history_note="file:"+filepath)
        return resultDict
    
    def addRowToMaster(self,row_to_add,history_note=""):
        self.master_history.append(str(len(self.master_history))+':'+history_note)
        self.master_tab.add_row(row_to_add)