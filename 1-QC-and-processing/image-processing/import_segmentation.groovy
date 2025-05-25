// Get image name
def entry = getProjectEntry()
def name = entry.getImageName()
name2 = name.substring(0,6)

def fileName = name2 + ".geojson"

// Build filepath
def pathInput = buildFilePath("/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/hcc-tma-project/baysor/qupath_project_post_baysor/hcc_tma2/geojsons", fileName)

// Import annotations
importObjectsFromFile(pathInput)