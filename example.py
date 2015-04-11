import numpy as np

from ota33a import OTA33A

exampleData = np.load('exampleData.npz')
gasConc = exampleData['gasConc']
temp = exampleData['temp']
pres = exampleData['pres']
ws3z = exampleData['ws3z']
ws3x = exampleData['ws3x']
ws3y = exampleData['ws3y']
distance = exampleData['distance']

a = OTA33A.fieldData(gasConc,temp,pres,ws3z,ws3x,ws3y,distance,16.04,'CH4')

massRate, volumeRate = a.getEmissionRate(make_plot=True)

print("The EPA's OTA33A method predicts an emission rate of %.2f g/s" % (massRate))
