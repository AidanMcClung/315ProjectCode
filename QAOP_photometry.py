from photutils.aperture import aperture_photometry #main photometry package/method - inherently linked to the apertures we import later
from photutils.aperture import ApertureStats #For getting the median and area that we use to calculate the background
from astropy.table import Table #For handling tables
from photutils.aperture import SkyCircularAperture, SkyCircularAnnulus
from astropy.coordinates import SkyCoord #for defining aperture positions
import astropy.units as u #Need 'deg' and 'arcsec' for skycoord and aperture
from astropy.io import fits #for loading fits files
from astropy.wcs import WCS #for getting WCS data from header
import os #for file handling and saving!
import io #for writing simple log files!

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

class photInstance:
    '''A class designed to be created in an external notebook and allow the easy use of the functionality of this module. When being created it will need to have a file of the apertures. By default it will assume this file is in the same root directory, but a path may be specified by passing it as the "apertureFilePath" parameter. An alternate directory for the module to store results can also be specified by passing the parameter "resultDir", which defaults to a folder called "photometry" in the root directory. 
    
    The primary method of the class is the runForFile() method, which takes in the filepath of the fits image it is to do the photometry one. The method automatically saves the result to a new row in the internal "Master Results Table", which is written to a file and saved after each call. As well, the method will return the dictionary of results, including the added Time column, should there be a need to use it elsewhere, though in general it is expected this will be discarded.
    
    As previously mentioned, the class has an internal "MasterTable" which stores the data in a 3D Table (Vertical is Time/File (it's meant to be time, but can have repeats), horizontal is sources (based on name), and lastly depth is the value stored (ie, aperture_raw_sum)). In the current implementation the three-dimensionality is not fully integrated into the table, and it is technically a 2D table where every cell is a dictionary. This has the drawback that one has to extract values cell-wise as it currently stands, but unpacking the data is handled by the "export" method and data can be more easily accesed there.
    
    The class also keeps a log of what it has run; it will write to this log "{IndexOfRowAdded}:{method}:{otherInformation}". As it currently stands the only possible method to add data is by using a fits file, so "file" is the only thing that will appear as "method". For files, the "otherInformation" is the path to this file. 
    (it is anticipated that future updates or external integration may use this log to avoid rerunning rows multiple times in the event of a runtime failure, but as of 07/28/2023 the info is just there and possible to be accesed by a user)
    
    NOTE: The class will open and load in previous results that are stored in the "master_table" and "master_log" files, as a way to continue in the event of a system failure. If it is desired that the entire process be started again, these files should be cleared/deleted manually before creating a new instance of this class.
    NOTE: This overwriting applies specifically and especially to any changes in the apertures. IF THE APERTURE FILE IS CHANGED between creations of the photInstance class, then the stored backup master_table file will have the WRONG NUMBER OF COLUMNS and the program will throw an error (or lots of them), so the backup must be cleared/deleted.'''
    #names = [] #I believe these C# style declarations are automatically handled by __init__.
    #apertures = []
    #annuli = []
    #master_history = []
    
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
            self.master_tab.write(self.resultDir+'/master_table.ecsv',format='ascii.ecsv')
            
        #okay, so I got too fancy, now I need to open up the log as well
        master_log_count = os.listdir(resultDir).count('master_log.txt')
        if master_log_count:
            with open(resultDir+'/master_log.txt','r') as log:
                self.master_history = log.readlines()
        #and lastly, for if/when we need them again, we add the parameters as variables to the class
        self.apertureFilePath,self.resultDir = apertureFilePath, resultDir
        #END INIT: Created self variables are [names,apertures,annuli,master_tab,master_history,apertureFilePath,resultDir]
    
    def createMasterTable(self):
        '''Uses the aperture names stored as "names" in the class instance to generate the master 3D table. Each row is an iteration of data addition (ie, a file that was read in and had photometry done on it), and each column is either an aperture, or, in the case of the first row, the datetime string for when the image was taken. The third dimension is acheived by storing a dictionary of the photometric results (raw sum, aperture area, local median, calculated local background, and adjusted sum) to each cell.
        
        This method can be called from anywhere to return an empty table with a column for each named source.
        
        NOTE: due to the implementation of how rows are added and such, changes to the apertures file after the photometric process has begun can cause errors, so any created files for the photometry process should be deleted if the apertures change.
        '''
        col_names = np.hstack(('Time',self.names)) #create a list of columns: Time | name1 | name2 | ...
        col_dtypes = np.hstack((str,np.full(len(self.names),dict))) #set first column to a string and the rest to dict: str | dict | dict | ...
        MasterResultTab = Table(names=col_names,dtype=col_dtypes)
        return MasterResultTab
    
    def runForFile(self,filepath):
        '''The primary method for the class, and the intended connection point between users and the module. It takes in a filepath as a parameter, and then does its photometry work. The method returns the dictionary it uses to add a row to the master table. It is expected that this is generally discarded, but it is provided should an external user ever find a need for it. 
        The method calls the modules doForApertures() method using the aperture lists it has stored internally, after extracting the datetime of the observation and WCS info, as well as the image, from the file it was provided.
        
        The method will save the results as a new row in its internal master table and add an entry in the log file containing the filepath it used.
        
        NOTE: This method does not discriminate, and will re-add rows as many times as it is called, so make sure to clear/delete the backup if you would like to avoid repeats. (And/Or filter them out afterward)'''
        image,wcs,img_time = loadImageAndWCS(filepath)
        resultDict = doForApertures(image,self.names,self.apertures,self.annuli,wcs)
        resultDict['Time'] = img_time
        self.addRowToMaster(resultDict,history_note="file:"+filepath)
        return resultDict
    
    def addRowToMaster(self,row_to_add,history_note=""):
        '''Adds a row to the internal master dictionary and an entry in the log. It is expecting to get a dictionary to add as the new row, with a named key/value pair for each column. It is untested what happens if you dont have this, but the author does not expect its a good thing at all, and did not integrate workarouds or checks originally because in their ideal world there are never errors and this is always being called by an instance with constant name/aperture lists.'''
        self.master_history.append(str(len(self.master_tab))+':'+history_note) #append a new status message to the log list
        self.master_tab.add_row(row_to_add) #add the row! it should have the right columns and column names otherwise it complains and dies
        #and lastly save the updated table
        self.master_tab.write(self.resultDir+'/master_table.ecsv',format='ascii.ecsv',overwrite=True)#save the new table
        with open(self.resultDir+'/master_log.txt','a') as log: #open in "append" mode
            log.write('\r\n'+self.master_history[-1]) #write the last item (which we added 4 lines above) to the file
        