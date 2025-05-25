// Cell Segmentation and Metadata Generation

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

// Get image name
def entry = getProjectEntry()
def name = entry.getImageName()
name2 = name.substring(0,6)

// Set pixel size
setPixelSizeMicrons(0.119831, 0.119831)

// For progress
print name2 + " DONE"
