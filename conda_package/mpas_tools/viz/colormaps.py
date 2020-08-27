import xml.etree.ElementTree as ET
import pkg_resources
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt


def register_sci_viz_colormaps():
    """Register all SciVisColor colormaps with matplotlib"""

    for mapName in ['3wave-yellow-grey-blue', '3Wbgy5',
                    '4wave-grey-red-green-mgreen', '5wave-yellow-brown-blue',
                    'blue-1', 'blue-3', 'blue-6', 'blue-8', 'blue-orange-div',
                    'brown-2', 'brown-5', 'brown-8', 'green-1', 'green-4',
                    'green-7', 'green-8', 'orange-5', 'orange-6',
                    'orange-green-blue-gray', 'purple-7', 'purple-8', 'red-1',
                    'red-3', 'red-4', 'yellow-1', 'yellow-7']:

        xmlFile = pkg_resources.resource_filename(
            __name__, 'SciVisColorColormaps/{}.xml'.format(mapName))
        _read_xml_colormap(xmlFile, mapName)


def _read_xml_colormap(xmlFile, mapName):
    """Read in an XML colormap"""

    xml = ET.parse(xmlFile)

    root = xml.getroot()
    colormap = root.findall('ColorMap')
    if len(colormap) > 0:
        colormap = colormap[0]
        colorDict = {'red': [], 'green': [], 'blue': []}
        for point in colormap.findall('Point'):
            x = float(point.get('x'))
            color = [float(point.get('r')), float(point.get('g')),
                     float(point.get('b'))]
            colorDict['red'].append((x, color[0], color[0]))
            colorDict['green'].append((x, color[1], color[1]))
            colorDict['blue'].append((x, color[2], color[2]))
        cmap = LinearSegmentedColormap(mapName, colorDict, 256)

        _register_colormap_and_reverse(mapName, cmap)


def _register_colormap_and_reverse(mapName, cmap):
    plt.register_cmap(mapName, cmap)
    plt.register_cmap('{}_r'.format(mapName), cmap.reversed())
