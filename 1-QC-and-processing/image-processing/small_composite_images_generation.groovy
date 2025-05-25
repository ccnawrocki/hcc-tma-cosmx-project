// Small Composite Images Generation
// Run without saving!

// Need this library
import qupath.lib.gui.images.servers.RenderedImageServer

// Set Channel Information and Colors
setChannelNames(
    "PanCK",
    "CK8/18", 
    "Membrane", 
    "CD45", 
    "DAPI"
)

setChannelColors(
    getColorRGB(0, 255, 0),
    getColorRGB(255, 255, 0), 
    getColorRGB(0, 255, 255), 
    getColorRGB(255, 0, 0), 
    getColorRGB(0, 0, 255)
)

// Set display ranges
createFullImageAnnotation(true)
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons":0.119831,"region":"ROI","tileSizeMicrons":510.0,"channel1":true,"channel2":true,"channel3":false,"channel4":true,"channel5":true,"doMean":true,"doStdDev":true,"doMinMax":false,"doMedian":true,"doHaralick":false}')

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()
def annotations = hierarchy.getAnnotationObjects()

vals = annotations.get(0).getMeasurements().values()

DAPI_max = vals[9] + 5*vals[10]
DAPI_min = vals[11] //+ 1*vals[10]
PanCK_max = vals[0] + 15*vals[1]
PanCK_min = vals[2] + 0*vals[1]
CD45_max = vals[6] + 10*vals[7]
CD45_min = vals[8] + 0*vals[7]
CK8_18_max = vals[3] + 6*vals[4]
CK8_18_min = vals[5] + 0*vals[4]

setChannelDisplayRange("Membrane", 255, 255)
setChannelDisplayRange("CK8/18", CK8_18_min, CK8_18_max)
setChannelDisplayRange("DAPI", DAPI_min, DAPI_max)
setChannelDisplayRange("PanCK", PanCK_min, PanCK_max)
setChannelDisplayRange("CD45", CD45_min, CD45_max)

// Get image name
def entry = getProjectEntry()
def name = entry.getImageName()
name2 = name.substring(0,6)

// Export small rendered RGB for plotting
def server = new RenderedImageServer.Builder(imageData)
    .downsamples(2)
    .build()

writeImage(server, "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc_tma1/small_rgb_images/" + name2 + ".jpeg")

