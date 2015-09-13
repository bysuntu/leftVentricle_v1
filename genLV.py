import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os
import csv
import numpy as np
from scipy import misc
import scipy.interpolate as si
from scipy.interpolate import pchip
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def itpl(PTS,numPTS,degree,sID,eID):
    t=range(len(PTS))
    ipl_t=np.linspace(sID,eID,numPTS*(eID-sID)/(len(PTS)-1)+1)
    newX=pchip(t,[row[0] for row in PTS])(ipl_t)
    newY=pchip(t,[row[1] for row in PTS])(ipl_t)
    newZ=pchip(t,[row[2] for row in PTS])(ipl_t)

    if degree==1:
        newX=interp1d(t,[row[0] for row in PTS],kind='slinear')(ipl_t)
        newY=interp1d(t,[row[1] for row in PTS],kind='slinear')(ipl_t)
        newZ=interp1d(t,[row[2] for row in PTS],kind='slinear')(ipl_t)
    elif degree==2:
        newX=interp1d(t,[row[0] for row in PTS],kind='quadratic')(ipl_t)
        newY=interp1d(t,[row[1] for row in PTS],kind='quadratic')(ipl_t)
        newZ=interp1d(t,[row[2] for row in PTS],kind='quadratic')(ipl_t)
    elif degree==3:
        newX=interp1d(t,[row[0] for row in PTS],kind='cubic')(ipl_t)
        newY=interp1d(t,[row[1] for row in PTS],kind='cubic')(ipl_t)
        newZ=interp1d(t,[row[2] for row in PTS],kind='cubic')(ipl_t)
    elif degree==4:
        newX=pchip(t,[row[0] for row in PTS])(ipl_t)
        newY=pchip(t,[row[1] for row in PTS])(ipl_t)
        newZ=pchip(t,[row[2] for row in PTS])(ipl_t)
    else:
        newX=interp1d(t,[row[0] for row in PTS],kind='quadratic')(ipl_t)
        newY=interp1d(t,[row[1] for row in PTS],kind='quadratic')(ipl_t)
        newZ=interp1d(t,[row[2] for row in PTS],kind='quadratic')(ipl_t)

    outPTS=[]
    for i in range(0,len(newX)):
        outPTS.append([newX[i],newY[i],newZ[i]])

    return outPTS

def longContour(sP, colorPixels,fName,frameNum,nPTS,pixelResolution):
    if len(sP)==0 or len(colorPixels)==0:
        print 'longContour Error: {}'.format(fName)
        return []

    #Interpolate
    ss=pixelResolution
    ll=pixelResolution+2
    
    nSP=[]
    oSP=sP
    nSP.append(sP)
    tail=[]
    while 1:       
        for i in range(0, len(colorPixels)):
            x=colorPixels[i][0]
            y=colorPixels[i][1]
            z=colorPixels[i][2]
            oDis=((oSP[0]-x)**2+(oSP[1]-y)**2+(oSP[2]-z)**2)**0.5
            nDis=((nSP[-1][0]-x)**2+(nSP[-1][1]-y)**2+(nSP[-1][2]-z)**2)**0.5
            if oDis>=nDis and nDis>ss and nDis<ll: # Not working when there is a sharp turning
                oSP=nSP[-1]
                nSP.append([x,y,z])
                break

        if i==(len(colorPixels)-1) and len(nSP)>=2: #Last iteration
            
            DIS=((nSP[-2][0]-nSP[-1][0])**2+(nSP[-2][1]-nSP[-1][1])**2+(nSP[-2][2]-nSP[-1][2])**2)**0.5
            for ii in range (0, len(colorPixels)):
                x=colorPixels[ii][0]
                y=colorPixels[ii][1]
                z=colorPixels[ii][2]
                oDis=((nSP[-2][0]-x)**2+(nSP[-2][1]-y)**2+(nSP[-2][2]-z)**2)**0.5
                nDis=((nSP[-1][0]-x)**2+(nSP[-1][1]-y)**2+(nSP[-1][2]-z)**2)**0.5
                if oDis>=nDis and oDis>DIS and nDis<=ss:
                    tail.append([x,y,z])
            break

    if len(tail)>0:
        oDis=((nSP[-1][0]-tail[0][0])**2+(nSP[-1][1]-tail[0][1])**2+(nSP[-1][2]-tail[0][2])**2)**0.5
        for i in range(0, len(tail)):
            nDis=((nSP[-1][0]-tail[i][0])**2+(nSP[-1][1]-tail[i][1])**2+(nSP[-1][2]-tail[i][2])**2)**0.5
            if nDis>=oDis:
                oDis=nDis
                id=i
        nSP.pop()
        nSP.append([tail[id][0],tail[id][1],tail[id][2]])
    #SPLINE
    numPTS=nPTS+1
    points=np.array(nSP)
    if len(points)>3:
        degree=3
    else:
        degree=len(points)-1
    x=points[:,0]
    y=points[:,1]
    z=points[:,2] 
    t=range(len(x))

    ipl_t=np.linspace(0,len(x)-1,numPTS)
    x_tup=si.splrep(t,x,k=degree,per=0)
    x_list=list(x_tup)
    xl=x.tolist()
    x_list[1]=xl
    x_i=si.splev(ipl_t,x_list)

    ipl_t=np.linspace(0,len(y)-1,numPTS)
    y_tup=si.splrep(t,y,k=degree,per=0)
    y_list=list(y_tup)
    yl=y.tolist()
    y_list[1]=yl
    y_i=si.splev(ipl_t,y_list)

    ipl_t=np.linspace(0,len(z)-1,numPTS)
    z_tup=si.splrep(t,z,k=degree,per=0)
    z_list=list(z_tup)
    zl=z.tolist()
    z_list[1]=zl
    z_i=si.splev(ipl_t,z_list)

    xyz=[]
    for i in range(0, len(x_i)):
        xyz.append([x_i[i],y_i[i],z_i[i]])
    xyz=itpl(points,200,degree,0,len(points)-1)
    return xyz

def loadContour(imageName, regionalColor, info):
    # Color
    colorMatrix = []
    colorMatrix.append(["RED",[1,0,0]]) #Red         0
    colorMatrix.append(["LIME",[0,1,0]]) #Lime        1
    colorMatrix.append(["BLUE",[0,0,1]]) #Blue        2
    colorMatrix.append(["GREEN",[0,0.5,0]]) #Green     3
    colorMatrix.append(["YELLOW",[1,1,0]]) #Yellow      4
    colorMatrix.append(["CYAN",[0,1,1]]) #Cyan        5
    colorMatrix.append(["MAGENTA",[1,0,1]]) #Magenta     6
    colorMatrix.append(["NAVY",[0,0,0.5]]) #Navy      7
    colorMatrix.append(["TEAL",[0,0.5,0.5]]) #Teal    8
    colorMatrix.append(["PURPLE",[0.5,0,0.5]]) #Purple  9
    colorMatrix.append(["OLIVE",[0.5,0.5,0]]) #Olive   10
    colorMatrix.append(["MAROON",[0.5,0,0]]) #Maroon    11
              
    numColors = len(regionalColor)
    colorID=[]
    for i in range(0, numColors):
        color = regionalColor[i]
        for j in range(0, 12):
            if color==colorMatrix[j][0]:
                colorID.append(j)
    
    colorPixels=[]
    for i in range(0, numColors):
        colorPixels.append([])
        
    if os.path.exists(imageName):
        fig = misc.imread(imageName)
        [r, c, oo] = np.shape(fig)
        
        for i in range(0, r):
            for j in range(0, c):
                for k in range(0, numColors):
                    colorIndex=colorMatrix[colorID[k]][1]
                    if colorMatrix[colorID[k]][0] == 'RED' and np.absolute(fig[i][j][0]-fig[i][j][1]) > 50 and np.absolute(fig[i][j][0]-fig[i][j][2]) > 50:
                        colorPixels[k].append([i,j])
                    if colorMatrix[colorID[k]][0] != 'RED' and np.absolute(fig[i][j][0]-(255*colorIndex[0])) < 20 and np.absolute(fig[i][j][1]-(255*colorIndex[1])) < 20 and np.absolute(fig[i][j][2]-int(255*colorIndex[2])) < 20:
                        colorPixels[k].append([i,j])
        return colorPixels
    else:
        return []

def startPt(pixels): #Calculating the center of certain colored region
    if len(pixels)==0:
        return []
    center=np.mean(pixels,axis=0)
    return center

def digitize3D(colorPixel, info):
    if colorPixel == []:
        return []
    origin=info[0:3]
    columnDir=info[3:6]
    rowDir=info[6:9]
    rS=info[11]
    cS=info[12]
    
    rotatedColor=[]
    for i in range(0, len(colorPixel)):
        rotatedColor.append(np.multiply(colorPixel[i][0]*rS, rowDir)+np.multiply(colorPixel[i][1]*cS,columnDir)+origin)
    return rotatedColor


def extractLongContours(i, info):
    twochName = "ch2\\ch2%04d.png" % i
    threechName = "ch3\\ch3%04d.png" % i
    fourchName = "ch4\\ch4%04d.png" % i

    # TwoChamber
    colors = ['RED', 'BLUE', 'LIME', 'GREEN']
    ch2Pixels = loadContour(twochName, colors, info[2])
    ch2PtsVen = digitize3D(ch2Pixels[0], info[2])
    ch2Blue = digitize3D(ch2Pixels[1], info[2])
    ch2Lime = digitize3D(ch2Pixels[2], info[2])
    # ch2Green = digitize3D(ch2Pixels[3], info[2])

    ch2Left = longContour(startPt(ch2Blue),
                                ch2PtsVen, 'ch2B', i, 100, 8)
    print '----------ch2Left'
    ch2Right = longContour(startPt(ch2Lime),
                                 ch2PtsVen, 'ch2G', i, 100, 8)
    print '----------ch2Right'

    print 'Loaded Two Chamber'

    # ThreeChamber
    colors = ['RED', 'BLUE', 'LIME', 'MAROON',
              'CYAN', 'MAGENTA', 'GREEN', 'YELLOW', 'OLIVE']
    ch3Pixels = loadContour(threechName, colors, info[3])
    ch3PtsVen = digitize3D(ch3Pixels[0], info[3])
    ch3Blue = digitize3D(ch3Pixels[1], info[3])
    ch3Lime = digitize3D(ch3Pixels[2], info[3])
    ch3Maroon = digitize3D(ch3Pixels[3], info[3])
    ch3Cyan = digitize3D(ch3Pixels[4], info[3])
    ch3Yellow = digitize3D(ch3Pixels[7], info[3])
    ch3Olive = digitize3D(ch3Pixels[8], info[3])

    ch3Left = longContour(startPt(ch3Blue),
                                ch3PtsVen, 'ch3B', i, 100, 8)
    print '----------ch3Left'
    ch3Right = longContour(startPt(ch3Lime),
                                 ch3PtsVen, 'ch3G', i, 100, 8)
    print '----------ch3Right'
    apex = startPt(ch3Maroon)
    venAortaAtrium = startPt(ch3Cyan)

    anterior = longContour(startPt(ch3Cyan),
                                 ch3Yellow, 'Anterior', i, 24, 4)
    print '----------anterior'
    posterior = longContour(startPt(ch3Lime),
                                  ch3Olive, 'Posterior', i, 24, 4)
    print '----------posterior'

    print 'Loaded Three Chamber'

    # FourChamber
    colors = ['RED', 'BLUE', 'LIME', 'GREEN']
    ch4Pixels = loadContour(fourchName, colors, info[4])
    ch4PtsVen = digitize3D(ch4Pixels[0], info[4])
    ch4Blue = digitize3D(ch4Pixels[1], info[4])
    ch4Lime = digitize3D(ch4Pixels[2], info[4])
    ch4Left = longContour(startPt(ch4Blue),
                                ch4PtsVen, 'ch4B', i, 100, 8)
    print '----------ch4Left'
    ch4Right = longContour(startPt(ch4Lime),
                                 ch4PtsVen, 'ch4G', i, 100, 8)
    print '----------ch4Right'

    print 'Loaded Four Chamber'

    bones = [ch2Left, ch2Right, ch3Left, ch3Right,
             ch4Left, ch4Right, anterior, posterior]
    keyPts = [apex, venAortaAtrium]

    return [bones, keyPts]

# Load Information
# fileName=raw_input('Key the DICOM information (csv): ')
fileName = 'MRIdata.csv'
info = np.zeros((40, 13))
rowNum = 0
with open(fileName, 'rb') as information:
    reader = csv.reader(information)
    for r in reader:
        for c in range(0, 14):
            if rowNum > 0 and c > 0:
                info[rowNum-1][c-1] = float(r[c])
        rowNum += 1

# Information
print 'Long Axis'
frameNum = int(raw_input('Total number of frames: '))
bSliceID = int(raw_input('   The bottom Slice ID: '))
frameNum = 40
bSliceID = 10

# Axis direction and apex based on frame0 !!!!!!!!!!!!!!
ooLong = extractLongContours(0, info)
keyPts = ooLong[1]
apex = keyPts[0]
# Low Slice and vertical Axis
lowSliceOrigin = info[bSliceID+5][0:3]
cDir = info[bSliceID+5][3:6]
rDir = info[bSliceID+5][6:9]

verticalDir = (keyPts[1]-keyPts[0])
if np.dot(verticalDir, np.cross(rDir, cDir)) > 0:
    axisDir = np.cross(rDir, cDir)
else:
    axisDir = -1*np.cross(rDir, cDir)
