import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os
import csv
import numpy as np
import sys

import p2dr1
import tools
import buildSTLr1


def writePts(pixels, frameNum):
    # Save file
    fileName = 'Pixel_{}.asc'.format(frameNum)
    file = open(fileName, "a")
    for i in range(0, len(pixels)):
        file.write('{}  {}  {} \n'.format(
            pixels[i][0], pixels[i][1], pixels[i][2]))
    file.close()
    return []


def reGenerateAtrium(atrialData, atrialVolume, axisDir, esID, id):
    atrialMESH = atrialData[0]
    inletRing = atrialMESH[len(atrialMESH)-1]

    refPt = np.mean(inletRing, axis=0)
    
    # Load STL files and calVolume
    if os.name == 'posix':
        stlVenName = 'Simulation/STL/atrium_{}.stl'.format(id)
        stlPosName = 'Simulation/STL/posterior_{}.stl'.format(id)
        stlAntName = 'Simulation/STL/anterior_{}'.stl.format(id)
    else:
        stlVenName = 'Simulation\\STL\\atrium_{}.stl'.format(id)
        stlPosName = 'Simulation\\STL\\posterior_{}.stl'.format(id)
        stlAntName = 'Simulation\\STL\\anterior_{}.stl'.format(id)
    stlNames = [stlVenName, stlPosName, stlAntName]
    vol = 0
    if id <= esID:
        for i in range(0, 3):
            stlName = stlNames[i]
            surface = tools.loadFacet(stlName)
            for n in range(0, len(surface)):
                facet = np.array(surface[n])
                vol += np.absolute(tools.tetraV(facet-refPt))
        aTip = atrialData[1]
        aFTip = aTip[0]
        aSTip = aTip[1]
        numSeg = len(aTip[0])/2
        pTip = atrialData[2]
        pFTip = pTip[0][0:len(pTip[0])-numSeg]
        pSTip = pTip[1][0:len(pTip[1])-numSeg]
        # TotalNumber of segments
        tNumSeg = 4*numSeg
        aFTip = tools.reGenLine([], aFTip, [], tNumSeg)
        aSTip = tools.reGenLine([], aSTip, [], tNumSeg)
        pFTip = tools.reGenLine([], pFTip, [], tNumSeg)
        pSTip = tools.reGenLine([], pSTip, [], tNumSeg)
        # Seal volume
        for i in range(0, tNumSeg):
            # First half
            p1 = aFTip[i]
            p2 = pFTip[i]
            p3 = aFTip[i+1]
            p4 = pFTip[i+1]
            facet1 = np.array([p1, p2, p3])
            facet2 = np.array([p2, p4, p3])
            vol += np.absolute(tools.tetraV(facet1-refPt))
            vol += np.absolute(tools.tetraV(facet2-refPt))
            # Second half
            p1 = aSTip[i]
            p2 = pSTip[i]
            p3 = aSTip[i+1]
            p4 = pSTip[i+1]
            facet1 = np.array([p1, p2, p3])
            facet2 = np.array([p2, p4, p3])
            vol += np.absolute(tools.tetraV(facet1-refPt))
            vol += np.absolute(tools.tetraV(facet2-refPt))
    else:
        for i in range(0, 1):
            stlName = stlNames[i]
            surface = tools.loadFacet(stlName)
            for n in range(0, len(surface)):
                facet = np.array(surface[n])
                vol += np.absolute(tools.tetraV(facet-refPt))
        jointRing = atrialMESH[0]
        jointCenter = np.mean(jointRing, axis=0)
        for i in range(0, len(jointRing)-1):
            p1 = jointRing[i]
            p2 = jointRing[i+1]
            facet = np.array([p1, p2, jointCenter])
            vol += np.absolute(tools.tetraV(facet-refPt))
        
    if atrialVolume < vol and id <= esID:
        sys.exit('Increase the atrial Volume')
    else:
        inletArea = 0
        for i in range(0, len(inletRing)-1):
            oo1 = inletRing[i]-refPt
            oo2 = inletRing[i+1]-refPt
            area = np.linalg.norm(np.cross(oo1, oo2))*0.5
            inletArea += area
        increaseVol = atrialVolume-vol
        increaseHeight = increaseVol/inletArea
        if increaseHeight < 0:
            increaseHeight = 0
    print 'Atrial volume: {} and inlet area: {}'.format(vol, inletArea)
    print 'Increase height: {}'.format(increaseHeight)
    newInletRing = []
    newAtrialMESH = []
    for i in range(0, len(atrialMESH[0])):
        pt = np.array(inletRing[i])+increaseHeight*axisDir
        newInletRing.append(pt)
        newLine = []
        for j in range(0, len(atrialMESH)):
            oo = atrialMESH[j][i]
            newLine.append(oo)
        if increaseHeight > 0:
            newLine.append(pt)
        finalLine = tools.reGenLine([], newLine, [], len(atrialMESH)-1)
        newAtrialMESH.append(finalLine[::-1])
    # New Atrium
    # Build STL File
    if os.name == 'posix':
        stlName = 'Simulation/STL/newAtrium_{}.stl'.format(id)
    else:
        stlName = 'Simulation\\STL\\newAtrium_{}.stl'.format(id)
    stlAtrium = open(stlName, "w")
    stlAtrium.write('solid Exported from mainFunc.py \n')
    tools.quadMesh(newAtrialMESH, [], [], [], [], [], [], [], stlAtrium, 'FULL')
    stlAtrium.write('endsolid Exported from topology.py')
    stlAtrium.close()

    # Inlet#
    # ------#
    # STL File
    if os.name == 'posix':
        stlName = 'Simulation/STL/newInlet_{}.stl'.format(id)
    else:
        stlName = 'Simulation\\STL\\newInlet_{}.stl'.format(id)

    stlInlet = open(stlName, "w")
    stlInlet.write('solid Exported from longAxis.py \n')
    topEdge = newInletRing
    tCenter = np.mean(topEdge, axis=0)
    centerEdge = []
    for i in range(0, len(topEdge)):
        centerEdge.append(tCenter)

    tools.quadMesh([centerEdge, topEdge], [], [], 20, [], [], 0, 1, stlInlet, 'Plane')
    stlInlet.write('endsolid Exported from topology.py')
    stlInlet.close()
    
    figure = plt.figure()
    ax = figure.gca(projection='3d')
    ax.scatter(refPt[0], refPt[1], refPt[2], c='b', marker='o', s=50)
    ax.plot([row[0] for row in newInletRing],
            [row[1] for row in newInletRing],
            [row[2] for row in newInletRing], '-ro')
    ax.plot([row[0] for row in inletRing],
            [row[1] for row in inletRing],
            [row[2] for row in inletRing], '-bo')
    if id > esID:
        ax.scatter(jointCenter[0], jointCenter[1], jointCenter[2], c='b', marker='o', s=50)
        ax.plot([row[0] for row in jointRing],
                [row[1] for row in jointRing],
                [row[2] for row in jointRing], '-gx')
    else:
        ax.plot([row[0] for row in aFTip],
                [row[1] for row in aFTip],
                [row[2] for row in aFTip],
                '-rx')
        ax.plot([row[0] for row in aSTip],
                [row[1] for row in aSTip],
                [row[2] for row in aSTip],
                '-ro')
        ax.plot([row[0] for row in pFTip],
                [row[1] for row in pFTip],
                [row[2] for row in pFTip],
                '-gx')
        ax.plot([row[0] for row in pSTip],
                [row[1] for row in pSTip],
                [row[2] for row in pSTip],
                '-go')
        for i in range(0, len(aSTip)):
            fLine = [aFTip[i], pFTip[i]]
            sLine = [aSTip[i], pSTip[i]]
            ax.plot([row[0] for row in fLine],
                    [row[1] for row in fLine],
                    [row[2] for row in fLine], '-')
    
    # plt.show()
    return 0

    
def reCoverMitral(angleLen, keyPts, axisDir, info, ax):
    # Load points
    keyPts = np.asarray(keyPts)
    joint = keyPts[1]
    p2 = keyPts[3]  # ch3Right
    p4 = keyPts[5]  # ch2Left
    p6 = keyPts[7]  # ch4Left
    topCenter = np.mean([p4, p6], axis=0)

    ax.scatter(joint[0], joint[1], joint[2], c='r', marker='o', s=100)
    ax.scatter(p2[0], p2[1], p2[2], c='g', marker='o', s=100)
    ax.scatter(topCenter[0], topCenter[1], topCenter[2],
               c='k', marker='x', s=50)
    # Angle between axisDir and ch3
    ch3Normal = np.cross(info[3][3:6], info[3][6:9])/np.linalg.norm(np.cross(
        info[3][3:6], info[3][6:9]))
    ooCh3 = axisDir-np.dot(ch3Normal, axisDir)*ch3Normal
    axisDirCh3 = ooCh3/np.linalg.norm(ooCh3)
    ooRadialCh3 = np.cross(ch3Normal, axisDir)
    radialDirCh3 = ooRadialCh3/np.linalg.norm(ooRadialCh3)
    # Show axisDirCh3
    oo1 = topCenter+10.0*axisDir
    axisDirLine = [topCenter, oo1]
    oo1 = topCenter+10.0*axisDirCh3
    axisDirCh3Line = [topCenter, oo1]
    oo1 = topCenter+10.0*radialDirCh3
    radialDirCh3Line = [topCenter, oo1]
    ax.plot([row[0] for row in axisDirLine],
            [row[1] for row in axisDirLine],
            [row[2] for row in axisDirLine], '-r')
    ax.plot([row[0] for row in axisDirCh3Line],
            [row[1] for row in axisDirCh3Line],
            [row[2] for row in axisDirCh3Line], '-g')
    ax.plot([row[0] for row in radialDirCh3Line],
            [row[1] for row in radialDirCh3Line],
            [row[2] for row in radialDirCh3Line], '-r')
    # Angle between axis and ch3
    ooAngle = np.arccos(np.dot(ch3Normal, axisDir))
    angleShortCh3 = 0.5*3.1415926-ooAngle
    # Posterior
    # Angle between projected axisDir onto ch3 and posterior
    ooPAngle = angleLen[0][1]
    pAngle = np.arccos(np.cos(ooPAngle/180.0*3.1415926)/np.cos(angleShortCh3))
    if ooPAngle/180.0*3.1415926+angleShortCh3 > 3.1415926:
        pAngle = 2*3.1415926-pAngle
    print pAngle
    # Radial Direction for posterior
    oo0 = joint-p2
    if np.dot(oo0, radialDirCh3) > 0:
        posteriorRadialDir = radialDirCh3
    else:
        posteriorRadialDir = radialDirCh3*-1

    oo1 = p2+10*posteriorRadialDir
    posteriorRadialDirLine = [p2, oo1]
    ax.plot([row[0] for row in posteriorRadialDirLine],
            [row[1] for row in posteriorRadialDirLine],
            [row[2] for row in posteriorRadialDirLine], '-r')
    # Posterior
    pLen = angleLen[0][2]
    deltaX = np.sin(pAngle)*pLen*posteriorRadialDir
    deltaY = np.cos(pAngle)*pLen*axisDirCh3
    pt = p2+deltaX+deltaY
    ax.scatter(pt[0], pt[1], pt[2], c='r', marker='<', s=100)
    oo0 = [p2, pt]
    posteriorLine = tools.reGenLine([], oo0, [], 50)
    ax.plot([row[0] for row in posteriorLine],
            [row[1] for row in posteriorLine],
            [row[2] for row in posteriorLine],
            '-x')
    # Anterior
    # Angle between projected axisDir onto ch3 and anterior
    ooAAngle = angleLen[1][1]
    aAngle = np.arccos(np.cos(ooAAngle/180.0*3.1415926)/np.cos(angleShortCh3))
    if ooAAngle/180.0*3.1415926+angleShortCh3 > 3.1415926:
        aAngle = 2*3.1415926-aAngle
    print aAngle
    # Radial Direction for Anterior
    anteriorRadialDir = -1.0*posteriorRadialDir
    oo1 = joint+10*anteriorRadialDir
    anteriorRadialDirLine = [joint, oo1]
    ax.plot([row[0] for row in anteriorRadialDirLine],
            [row[1] for row in anteriorRadialDirLine],
            [row[2] for row in anteriorRadialDirLine], '-g')    
    # Anterior
    aLen = angleLen[1][2]
    deltaX = np.sin(aAngle)*aLen*anteriorRadialDir
    deltaY = np.cos(aAngle)*aLen*axisDirCh3
    pt = joint+deltaX+deltaY
    ax.scatter(pt[0], pt[1], pt[2], c='r', marker='>', s=100)
    oo0 = [joint, pt]
    anteriorLine = tools.reGenLine([], oo0, [], 50)
    ax.plot([row[0] for row in anteriorLine],
            [row[1] for row in anteriorLine],
            [row[2] for row in anteriorLine],
            '-+')
    return [posteriorLine, anteriorLine]


def projectPosterior(ch3Posterior, venLine, axisDir, topCenter, ax):
    # figure = plt.figure()
    # ax = figure.gca(projection='3d')
    # ch3 radial line
    oo0 = np.cross(axisDir, (ch3Posterior[0]-topCenter))
    ooDir = oo0/np.linalg.norm(oo0)
    oo0 = np.cross(axisDir, ooDir)
    ch3RadialDir = oo0/np.linalg.norm(oo0)
    # ven radial line
    ooDir = np.cross(axisDir, (venLine[0]-topCenter))/np.linalg.norm(np.cross(
        axisDir, (venLine[0]-topCenter)))
    radialDir = np.cross(axisDir, ooDir)/np.linalg.norm(np.cross(axisDir,
                                                                 ooDir))
    # radial line

    aList = [0]
    lList = [0]

    newPosterior = [venLine[0]]
    for i in range(1, len(ch3Posterior)):
        oo = np.array(ch3Posterior[i])-np.array(ch3Posterior[0])
        lenOO = np.linalg.norm(oo)
        lenAxisDir = np.linalg.norm(axisDir)
        ooAngle = np.arccos(np.dot(oo, axisDir)/(lenOO*lenAxisDir))
        if np.dot(oo, ch3RadialDir) >= 0:
            angle = ooAngle
        else:
            angle = 2.0*3.1415926-ooAngle
        length = np.linalg.norm(oo)
        aList.append(angle)
        lList.append(length)
        pt = length*np.cos(angle)*axisDir+length*np.sin(angle)*radialDir
        pt = pt+venLine[0]
        newPosterior.append(pt)

    # fakePosterior along venline
    fPosterior = [venLine[0]]
    for n in range(1, 5):
        seg = (venLine[n]-venLine[0])
        lenSeg = np.linalg.norm(seg)
        ooAngle = np.arccos(np.dot(seg, axisDir)/(lenSeg*lenAxisDir))
        if np.dot(seg, radialDir) >= 0:
            segAngle = ooAngle
        else:
            segAngle = 2.0*3.1415926-ooAngle
        shiftAngle = segAngle-1.0/3.0*3.1415926
        pt = lenSeg*np.cos(shiftAngle)*axisDir+lenSeg*np.sin(shiftAngle)*radialDir
        pt = pt+venLine[0]
        delDis = pt-(radialDir*lenSeg+venLine[n])
        fPosterior.append(pt)
    for i in range(5, len(ch3Posterior)):
        # lenSeg = np.linalg.norm(venLine[i]-venLine[i-1])
        pt = radialDir*lenSeg+venLine[i]+delDis
        fPosterior.append(pt)

    # Avoid touching
    checkDir = np.cross(axisDir, radialDir)
    unPosterior = [newPosterior[0]]
    for i in range(1, len(ch3Posterior)):
        ooLeafletLine = newPosterior[i]-newPosterior[0]
        ooVenLine = fPosterior[i]-fPosterior[0]
        ooDir = np.cross(ooLeafletLine, ooVenLine)
        if np.dot(checkDir, ooDir) < 0:
            oo0 = newPosterior[i:len(newPosterior)]-newPosterior[i]
            oo1 = oo0+fPosterior[i]
            newPosterior = tools.combineArrays([fPosterior[0:i], oo1])
            unPosterior.append(fPosterior[i])
        else:
            unPosterior.append(newPosterior[i])
    '''
    ax.plot([row[0] for row in newPosterior], [row[1] for row in newPosterior],
            [row[2] for row in newPosterior], '-r')
    ax.plot([row[0] for row in ch3Posterior], [row[1] for row in ch3Posterior],
            [row[2] for row in ch3Posterior], '-o')
    ax.plot([row[0] for row in fPosterior], [row[1] for row in fPosterior],
            [row[2] for row in fPosterior], '-go')
    
    ax.plot([row[0] for row in unPosterior], [row[1] for row in unPosterior],
            [row[2] for row in unPosterior], '-g')
    ax.plot([row[0] for row in venLine], [row[1] for row in venLine],
            [row[2] for row in venLine], '-r')
    '''
    '''
    rLine = [venLine[0], venLine[0]+10*radialDir]
    vLine = [venLine[0], venLine[0]+10*axisDir]
    chRLine = [ch3Posterior[0], ch3Posterior[0]+10*ch3RadialDir]
    ax.plot([row[0] for row in rLine], [row[1] for row in rLine],
            [row[2] for row in rLine], '-r')
    ax.plot([row[0] for row in vLine], [row[1] for row in vLine],
            [row[2] for row in vLine], '-b')
    ax.plot([row[0] for row in chRLine], [row[1] for row in chRLine],
            [row[2] for row in chRLine], '-o')
    '''
    # plt.show()
    return unPosterior


def GEO(ratio, keyPts, Bones, mitral, info, axisDir, esID, name, ax):
    # Ventricle
    venBone = Bones[0]
    # Points
    cutPt = (np.array(keyPts[1])-np.array(keyPts[0]))*ratio+np.array(keyPts[0])
    apex = np.array(keyPts[0])
    joint = np.array(keyPts[1])
    # Long contours
    ch2Left = venBone[0]
    ch2Right = venBone[1]
    ch3Left = venBone[2]
    ch3Right = venBone[3]
    ch4Left = venBone[4]
    ch4Right = venBone[5]
    l3L2L = venBone[6]
    l2L4R = venBone[7]
    l4R3R = venBone[8]
    l3R2R = venBone[9]
    l2R4L = venBone[10]
    l4L3L = venBone[11]

    # Upper ventricle
    # numSection=16
    numSeg = 6
    # --------------#
    # TopSlice and Ventricle Edge
    topCenter = (np.array(ch2Left[0])+np.array(ch4Left[0]))*0.5
    topSliceP1 = tools.l2s(ch2Right, cutPt, axisDir)
    topSliceP2 = tools.l2s(ch3Right, cutPt, axisDir)
    topSliceP3 = tools.l2s(ch4Right, cutPt, axisDir)
    topSliceP4 = tools.l2s(ch2Left, cutPt, axisDir)
    topSliceP5 = tools.l2s(ch3Left, cutPt, axisDir)
    topSliceP6 = tools.l2s(ch4Left, cutPt, axisDir)

    topSliceC = np.mean([topSliceP1, topSliceP2, topSliceP3,
                         topSliceP4, topSliceP5, topSliceP6], axis=0)
    # venZAxis=(topCenter-topSliceC)/
    # np.linalg.norm(topCenter-topSliceC)
    # venYAxis=np.cross(venZAxis,np.array(ch3Left[0])-np.array(topCenter))/
    # np.linalg.norm(np.cross(venZAxis,np.array(ch3Left[0])-np.array(topCenter)))
    # venXAxis=np.cross(venYAxis,venZAxis)/np.linalg.norm(np.cross(venYAxis,venZAxis))

    # Split the ventricle
    lowPts = 20
    highPts = 10
    ch3Right.reverse()
    ch4Right.reverse()
    ch2Left.reverse()
    ch3Left.reverse()
    ch4Left.reverse()
    ch2Right.reverse()
    l3L2L.reverse()
    l2L4R.reverse()
    l4R3R.reverse()
    l3R2R.reverse()
    l2R4L.reverse()
    l4L3L.reverse()

    lowCh2Left = tools.trimLine(ch2Left, topSliceC, axisDir, lowPts)
    lowCh2Right = tools.trimLine(ch2Right, topSliceC, axisDir, lowPts)
    lowCh3Left = tools.trimLine(ch3Left, topSliceC, axisDir, lowPts)
    lowCh3Right = tools.trimLine(ch3Right, topSliceC, axisDir, lowPts)
    lowCh4Left = tools.trimLine(ch4Left, topSliceC, axisDir, lowPts)
    lowCh4Right = tools.trimLine(ch4Right, topSliceC, axisDir, lowPts)
    low3L2L = tools.trimLine(l3L2L, topSliceC, axisDir, lowPts)
    low2L4R = tools.trimLine(l2L4R, topSliceC, axisDir, lowPts)
    low4R3R = tools.trimLine(l4R3R, topSliceC, axisDir, lowPts)
    low3R2R = tools.trimLine(l3R2R, topSliceC, axisDir, lowPts)
    low2R4L = tools.trimLine(l2R4L, topSliceC, axisDir, lowPts)
    low4L3L = tools.trimLine(l4L3L, topSliceC, axisDir, lowPts)

    lowCh3Right.reverse()
    lowCh4Right.reverse()
    lowCh2Left.reverse()
    lowCh3Left.reverse()
    lowCh4Left.reverse()
    lowCh2Right.reverse()
    low3L2L.reverse()
    low2L4R.reverse()
    low4R3R.reverse()
    low3R2R.reverse()
    low2R4L.reverse()
    low4L3L.reverse()
    lowVenMESH = [lowCh3Right, low4R3R, lowCh4Right, low2L4R, lowCh2Left,
                  low3L2L, lowCh3Left, low4L3L, lowCh4Left, low2R4L,
                  lowCh2Right, low3R2R]

    ch3Right.reverse()
    ch4Right.reverse()
    ch2Left.reverse()
    ch3Left.reverse()
    ch4Left.reverse()
    ch2Right.reverse()
    l3L2L.reverse()
    l2L4R.reverse()
    l4R3R.reverse()
    l3R2R.reverse()
    l2R4L.reverse()
    l4L3L.reverse()

    ch2Left = tools.trimLine(ch2Left, topSliceC, axisDir, highPts)
    ch2Right = tools.trimLine(ch2Right, topSliceC, axisDir, highPts)
    ch3Left = tools.trimLine(ch3Left, topSliceC, axisDir, highPts)
    ch3Right = tools.trimLine(ch3Right, topSliceC, axisDir, highPts)
    ch4Left = tools.trimLine(ch4Left, topSliceC, axisDir, highPts)
    ch4Right = tools.trimLine(ch4Right, topSliceC, axisDir, highPts)
    l3L2L = tools.trimLine(l3L2L, topSliceC, axisDir, highPts)
    l2L4R = tools.trimLine(l2L4R, topSliceC, axisDir, highPts)
    l4R3R = tools.trimLine(l4R3R, topSliceC, axisDir, highPts)
    l3R2R = tools.trimLine(l3R2R, topSliceC, axisDir, highPts)
    l2R4L = tools.trimLine(l2R4L, topSliceC, axisDir, highPts)
    l4L3L = tools.trimLine(l4L3L, topSliceC, axisDir, highPts)

    ringPTS = [ch3Right[0], l4R3R[0], ch4Right[0], l2L4R[0], ch2Left[0],
               l3L2L[0], ch3Left[0], l4L3L[0], ch4Left[0], l2R4L[0],
               ch2Right[0], l3R2R[0]]
    venRing = tools.shortAxis(ringPTS, topCenter, 1, len(ringPTS)*numSeg,
                              0, len(ringPTS), axisDir)
    # Normal Angle
    venNormalAngle = []
    for i in range(0, len(venRing)):
        oo = np.array(venRing[i])-np.array(topCenter)
        venNormalAngle.append(np.cross(axisDir, oo)/np.linalg.norm(
            np.cross(axisDir, oo)))

    # Building
    upVenMESH = []
    # Mitral leaflets
    figure = plt.figure()
    ax = figure.gca(projection='3d')
    ooPosterior = mitral[0]
    ooAnterior = mitral[1]
    ch3Posterior = tools.reGenLine([], ooPosterior, [], 4*numSeg)
    ch3Anterior = tools.reGenLine([], ooAnterior, [], 4*numSeg)
    posteriorMESH = []

    guidingRing = []
    for n in range(0, len(venNormalAngle)):
        ooMesh = [venRing[n]]
        for i in range(1, len(ch3Right)):
            p1 = ch2Right[i]
            p2 = ch3Right[i]
            p3 = ch4Right[i]
            p4 = ch2Left[i]
            p5 = ch3Left[i]
            p6 = ch4Left[i]
            l12 = l3R2R[i]
            l23 = l4R3R[i]
            l34 = l2L4R[i]
            l45 = l3L2L[i]
            l56 = l4L3L[i]
            l61 = l2R4L[i]
            venPTS = [p2, l23, p3, l34, p4, l45, p5, l56, p6, l61, p1, l12]
            ooRing = tools.shortAxis(venPTS, topCenter, 1, len(venPTS)*numSeg,
                                     0, len(venPTS), axisDir)
            oo = tools.l2sPt(ooRing, topCenter, venNormalAngle[n], axisDir)
            ooMesh.append(oo)
        if ch3Posterior != []:
            ooPVen = tools.reGenLine(ch3Posterior, ooMesh, [], [])
            if len(ooPVen) != len(ch3Posterior):
                print len(ooPVen)
                print len(ch3Posterior)
                '''
                ax.plot([row[0] for row in ooPVen],
                        [row[1] for row in ooPVen],
                        [row[2] for row in ooPVen], '-x')
                ax.plot([row[0] for row in ch3Posterior],
                        [row[1] for row in ch3Posterior],
                        [row[2] for row in ch3Posterior], '-o')
                '''
                sys.exit('The ratio is to big')
            newPosteriorLine = projectPosterior(ch3Posterior,
                                                ooPVen, axisDir, topCenter, ax)
            posteriorMESH.append(newPosteriorLine)
            '''
            ax.plot([row[0] for row in newPosteriorLine],
                    [row[1] for row in newPosteriorLine],
                    [row[2] for row in newPosteriorLine], '-kx')
            '''
        '''
        ax.plot([row[0] for row in ooMesh], [row[1] for row in ooMesh],
                [row[2] for row in ooMesh], '-r')
        '''
        upVenMESH.append(ooMesh)
        guidingRing.append(oo)

    buildSTLr1.ventricle(upVenMESH, lowVenMESH, topSliceC, apex, joint,
                         guidingRing, axisDir, len(venPTS), numSeg, name, ax)

    ax.plot([row[0] for row in lowCh3Right], [row[1] for row in lowCh3Right],
            [row[2] for row in lowCh3Right], '-x')
    ax.plot([row[0] for row in lowCh3Left], [row[1] for row in lowCh3Left],
            [row[2] for row in lowCh3Left], '-x')
    ax.plot([row[0] for row in lowCh2Left], [row[1] for row in lowCh2Left],
            [row[2] for row in lowCh2Left], '-x')
    ax.plot([row[0] for row in lowCh2Right], [row[1] for row in lowCh2Right],
            [row[2] for row in lowCh2Right], '-x')
    ax.plot([row[0] for row in lowCh4Left], [row[1] for row in lowCh4Left],
            [row[2] for row in lowCh4Left], '-x')
    ax.plot([row[0] for row in lowCh4Right], [row[1] for row in lowCh4Right],
            [row[2] for row in lowCh4Right], '-x')

    ax.plot([row[0] for row in ch3Right], [row[1] for row in ch3Right],
            [row[2] for row in ch3Right], '-o')
    ax.plot([row[0] for row in ch3Left], [row[1] for row in ch3Left],
            [row[2] for row in ch3Left], '-o')
    ax.plot([row[0] for row in ch2Left], [row[1] for row in ch2Left],
            [row[2] for row in ch2Left], '-o')
    ax.plot([row[0] for row in ch2Right], [row[1] for row in ch2Right],
            [row[2] for row in ch2Right], '-o')
    ax.plot([row[0] for row in ch4Left], [row[1] for row in ch4Left],
            [row[2] for row in ch4Left], '-o')
    ax.plot([row[0] for row in ch4Right], [row[1] for row in ch4Right],
            [row[2] for row in ch4Right], '-o')

    # Atrium
    atrialBone = Bones[1]
    atriumFull = buildSTLr1.atrium(atrialBone, venRing,
                                   axisDir, numSeg, name, ax)
    mitralRing = atriumFull[0]
    atrialMESH = atriumFull[1]
    # Mitral
    # Posterior

    if ch3Posterior != []:
        posteriorTip = buildSTLr1.posterior(posteriorMESH, topCenter,
                                            numSeg, axisDir,
                                            name, ax)
        scaledPosterior = posteriorTip[4]

        firstHalf = posteriorTip[0]
        secondHalf = posteriorTip[1]
        firstTop = posteriorTip[2]
        secondTop = posteriorTip[3]
        posteriorTop = [firstTop, secondTop]
        ax.plot([row[0] for row in firstHalf], [row[1] for row in firstHalf],
                [row[2] for row in firstHalf], '-x')
        ax.plot([row[0] for row in secondHalf], [row[1] for row in secondHalf],
                [row[2] for row in secondHalf], '-o')
        ax.plot([row[0] for row in firstTop],
                [row[1] for row in firstTop],
                [row[2] for row in firstTop], '-')
        ax.plot([row[0] for row in secondTop],
                [row[1] for row in secondTop],
                [row[2] for row in secondTop], '-')
        ax.scatter(secondHalf[0][0], secondHalf[0][1],
                   secondHalf[0][2], c='r', marker='>', s=100)
        ax.scatter(firstHalf[0][0], firstHalf[0][1],
                   firstHalf[0][2], c='r', marker='<', s=100)

    # Anterior
    if ch3Anterior != [] and ch3Posterior != []:
        aTip = buildSTLr1.anterior(ch3Anterior, scaledPosterior, posteriorTip, posteriorTop,
                                   mitralRing, axisDir, numSeg, esID, name, ax)

    # Aorta
    aorticBone = Bones[2]
    aorticZAxis = Bones[3]
    buildSTLr1.aorta(aorticBone, venRing, mitralRing,
                     aorticZAxis, numSeg, name, ax)

    return [atrialMESH, aTip, [firstHalf, secondHalf]]


def restoreGeo(axisDir, keyPts, centerLine, shorts, bSliceID, info, ax):
    keyPts = np.asarray(keyPts)
    apex = keyPts[0]
    joint = keyPts[1]
    p1 = keyPts[2]  # ch2Right
    p2 = keyPts[3]  # ch3Right
    p3 = keyPts[4]  # ch4Right
    p4 = keyPts[5]  # ch2Left
    p5 = keyPts[6]  # ch3Left
    p6 = keyPts[7]  # ch4Left
    venPTS = [p2, p3, p4, p5, p6, p1]
    topCenter = np.mean([p4, p6], axis=0)

    venRing = tools.shortAxis(venPTS, topCenter, 1, 120, 0,
                              len(venPTS), axisDir)
    shortCenter = []
    for i in range(0, bSliceID+1):
        if shorts[0][i] != []:
            center = np.mean(shorts[0][i], axis=0)
            ax.scatter(center[0], center[1], center[2], c='b', marker='x', s=50)
        else:
            center = []
        shortCenter.append(center)
    # Restore the short axis contours based on the long axis center line
    # Project the long axis contours onto the short axis
    # Contours
    nCh2Left = [p4]
    nCh2Right = [p1]
    nCh3Left = [p5]
    nCh3Right = [p2]
    nCh4Left = [p6]
    nCh4Right = [p3]
    # Additional lines
    ooRing = tools.shortAxis(venPTS, topCenter, 1, 12, 0, len(venPTS), axisDir)
    l3L2L = [ooRing[5]]
    l2L4R = [ooRing[3]]
    l4R3R = [ooRing[1]]
    l3R2R = [ooRing[11]]
    l2R4L = [ooRing[9]]
    l4L3L = [ooRing[7]]

    # Normal
    ch2Normal = np.cross(info[2][3:6], info[2][6:9])/np.linalg.norm(np.cross(
        info[2][3:6], info[2][6:9]))
    ch3Normal = np.cross(info[3][3:6], info[3][6:9])/np.linalg.norm(np.cross(
        info[3][3:6], info[3][6:9]))
    ch4Normal = np.cross(info[4][3:6], info[4][6:9])/np.linalg.norm(np.cross(
        info[4][3:6], info[4][6:9]))
    # Origin
    ch2Origin = np.mean([p4, p1], axis=0)
    ch3Origin = np.mean([p5, p2], axis=0)
    ch4Origin = np.mean([p6, p3], axis=0)
    nShorts = [venRing]
    for i in range(0, bSliceID+1):
        refPt = info[i+5][0:3]
        check1 = (np.dot(p1-refPt, axisDir) >= 0)
        check2 = (np.dot(p2-refPt, axisDir) >= 0)
        check3 = (np.dot(p3-refPt, axisDir) >= 0)
        check4 = (np.dot(p4-refPt, axisDir) >= 0)
        check5 = (np.dot(p5-refPt, axisDir) >= 0)
        check6 = (np.dot(p6-refPt, axisDir) >= 0)
        if check1 and check2 and check3 and check4 and check5 and check6:
            check = 1
        else:
            check = 0

        longCenterPt = tools.l2s(centerLine, refPt, axisDir)
        # Direction !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        flip = -1.0
        if longCenterPt != [] and shorts[0][i] != [] and check == 1:
            '''
            ax.scatter(longCenterPt[0], longCenterPt[1],
                       longCenterPt[2], c='r',
                       marker='v', s=100)
            '''
            dif = shortCenter[i]-longCenterPt
            ooShort = np.array(shorts[0][i])-np.array(dif)
            nShorts.append(ooShort)
            ax.plot([row[0] for row in ooShort], [row[1] for row in ooShort],
                    [row[2] for row in ooShort], '-r')
            # PAY ATTENTION ORIENTATION !!!!!!!!!!!!
            o2Left = tools.l2sPt(ooShort, ch2Origin, ch2Normal*flip, axisDir)
            nCh2Left.append(o2Left)
            o2Right = tools.l2sPt(ooShort, ch2Origin, -1*ch2Normal*flip, axisDir)
            nCh2Right.append(o2Right)
            o3Left = tools.l2sPt(ooShort, ch3Origin, ch3Normal*flip, axisDir)
            nCh3Left.append(o3Left)
            o3Right = tools.l2sPt(ooShort, ch3Origin, -1*ch3Normal*flip, axisDir)
            nCh3Right.append(o3Right)
            o4Left = tools.l2sPt(ooShort, ch4Origin, -1*ch4Normal*flip, axisDir)
            nCh4Left.append(o4Left)
            o4Right = tools.l2sPt(ooShort, ch4Origin, ch4Normal*flip, axisDir)
            nCh4Right.append(o4Right)
            # Additional
            l3L2L.append(tools.midPt(ooShort, o3Left, o2Left))
            l2L4R.append(tools.midPt(ooShort, o2Left, o4Right))
            l4R3R.append(tools.midPt(ooShort, o4Right, o3Right))
            l3R2R.append(tools.midPt(ooShort, o3Right, o2Right))
            l2R4L.append(tools.midPt(ooShort, o2Right, o4Left))
            l4L3L.append(tools.midPt(ooShort, o4Left, o3Left))

    ax.plot([row[0] for row in nCh2Left], [row[1] for row in nCh2Left],
            [row[2] for row in nCh2Left], '-r')
    ax.plot([row[0] for row in nCh2Right], [row[1] for row in nCh2Right],
            [row[2] for row in nCh2Right], '-xr')
    ax.plot([row[0] for row in nCh3Left], [row[1] for row in nCh3Left],
            [row[2] for row in nCh3Left], '-g')
    ax.plot([row[0] for row in nCh3Right], [row[1] for row in nCh3Right],
            [row[2] for row in nCh3Right], '-xg')
    ax.plot([row[0] for row in nCh4Left], [row[1] for row in nCh4Left],
            [row[2] for row in nCh4Left], '-b')
    ax.plot([row[0] for row in nCh4Right], [row[1] for row in nCh4Right],
            [row[2] for row in nCh4Right], '-xb')
    # plt.show()
    ax.plot([row[0] for row in l3L2L], [row[1] for row in l3L2L],
            [row[2] for row in l3L2L], '-')
    ax.plot([row[0] for row in l2L4R], [row[1] for row in l2L4R],
            [row[2] for row in l2L4R], '-')
    ax.plot([row[0] for row in l4R3R], [row[1] for row in l4R3R],
            [row[2] for row in l4R3R], '-')
    ax.plot([row[0] for row in l3R2R], [row[1] for row in l3R2R],
            [row[2] for row in l3R2R], '-')
    ax.plot([row[0] for row in l2R4L], [row[1] for row in l2R4L],
            [row[2] for row in l2R4L], '-')
    ax.plot([row[0] for row in l4L3L], [row[1] for row in l4L3L],
            [row[2] for row in l4L3L], '-')
    
    numPTS = 20
    # BSpline Interpolation
    iCh2Left = tools.bSLine(nCh2Left, 4, numPTS)
    iCh2Right = tools.bSLine(nCh2Right, 4, numPTS)
    iCh3Left = tools.bSLine(nCh3Left, 4, numPTS)
    iCh3Right = tools.bSLine(nCh3Right, 4, numPTS)
    iCh4Left = tools.bSLine(nCh4Left, 4, numPTS)
    iCh4Right = tools.bSLine(nCh4Right, 4, numPTS)

    ax.plot([row[0] for row in iCh2Left], [row[1] for row in iCh2Left],
            [row[2] for row in iCh2Left], '-or')
    ax.plot([row[0] for row in iCh2Right], [row[1] for row in iCh2Right],
            [row[2] for row in iCh2Right], '-vr')
    ax.plot([row[0] for row in iCh3Left], [row[1] for row in iCh3Left],
            [row[2] for row in iCh3Left], '-og')
    ax.plot([row[0] for row in iCh3Right], [row[1] for row in iCh3Right],
            [row[2] for row in iCh3Right], '-vg')
    ax.plot([row[0] for row in iCh4Left], [row[1] for row in iCh4Left],
            [row[2] for row in iCh4Left], '-ob')
    ax.plot([row[0] for row in iCh4Right], [row[1] for row in iCh4Right],
            [row[2] for row in iCh4Right], '-vb')

    iCh3L2L = tools.bSLine(l3L2L, 4, numPTS)
    iCh2L4R = tools.bSLine(l2L4R, 4, numPTS)
    iCh4R3R = tools.bSLine(l4R3R, 4, numPTS)
    iCh3R2R = tools.bSLine(l3R2R, 4, numPTS)
    iCh2R4L = tools.bSLine(l2R4L, 4, numPTS)
    iCh4L3L = tools.bSLine(l4L3L, 4, numPTS)

    ax.plot([row[0] for row in iCh3L2L], [row[1] for row in iCh3L2L],
            [row[2] for row in iCh3L2L], '-o')
    ax.plot([row[0] for row in iCh2L4R], [row[1] for row in iCh2L4R],
            [row[2] for row in iCh2L4R], '-v')
    ax.plot([row[0] for row in iCh4R3R], [row[1] for row in iCh4R3R],
            [row[2] for row in iCh4R3R], '-o')
    ax.plot([row[0] for row in iCh3R2R], [row[1] for row in iCh3R2R],
            [row[2] for row in iCh3R2R], '-v')
    ax.plot([row[0] for row in iCh2R4L], [row[1] for row in iCh2R4L],
            [row[2] for row in iCh2R4L], '-o')
    ax.plot([row[0] for row in iCh4L3L], [row[1] for row in iCh4L3L],
            [row[2] for row in iCh4L3L], '-v')

    # Atrium
    ooPts = [iCh2Left[0], joint, iCh4Left[0],
             iCh2Right[0], iCh3Right[0], iCh4Right[0]]
    ooCenter = np.mean(ooPts, axis=0)
    mitralRing = tools.shortAxis(ooPts, ooCenter, 1, len(ooPts)*10,
                                 0, len(ooPts), axisDir)

    ax.plot([row[0] for row in mitralRing], [row[1] for row in mitralRing],
            [row[2] for row in mitralRing], '-g')

    avPts = [iCh3Right[0], iCh4Right[0], iCh2Left[0],
             joint, iCh4Left[0], iCh2Right[0]]

    atrialPts = []
    atrialShorts = []
    for i in reversed(range(0, len(shorts[0]))):
        atrialShort = shorts[0][i]
        if atrialShort != []:
            center = np.mean(atrialShort, axis=0)
            sliceTop2Apex = tools.l2s([apex, joint], center, axisDir)
            if sliceTop2Apex == []:
                ax.plot([row[0] for row in atrialShort],
                        [row[1] for row in atrialShort],
                        [row[2] for row in atrialShort], '-vg')
                ax.scatter(center[0], center[1], center[2],
                           c='r', marker='o', s=50)

                ooPts = []
                for j in range(0, len(avPts)):
                    ooPts.append(tools.p2L(avPts[j], atrialShort))
                atrialPts.append(ooPts)
                atrialShorts.append(atrialShort)
    if atrialPts == []:
        sys.exit('No atrium')
    atrial3Right = [iCh3Right[0]]
    atrial4R3R = [iCh4R3R[0]]  #
    atrial4Right = [iCh4Right[0]]
    atrial2L4R = [iCh2L4R[0]]  #
    atrial2Left = [iCh2Left[0]]
    atrial2LJ = [tools.midPt(mitralRing, iCh2Left[0], joint)]  #
    atrialJoint = [joint]
    atrialJ4L = [tools.midPt(mitralRing, joint, iCh4Left[0])]  #
    atrial4Left = [iCh4Left[0]]
    atrial2R4L = [iCh2R4L[0]]  #
    atrial2Right = [iCh2Right[0]]
    atrial3R2R = [iCh3R2R[0]]  #

    for i in range(0, len(atrialPts)):
        atrialShort = atrialShorts[i]
        # atrial3Right.append(atrialPts[i][0])
        # atrial4Right.append(atrialPts[i][1])
        # atrial2Left.append(atrialPts[i][2])
        # atrialJoint.append(atrialPts[i][3])
        # atrial4Left.append(atrialPts[i][4])
        # atrial2Right.append(atrialPts[i][5])
        atrial4R3R.append(tools.midPt(atrialShort,
                                      atrialPts[i][0], atrialPts[i][1]))
        atrial2L4R.append(tools.midPt(atrialShort,
                                      atrialPts[i][1], atrialPts[i][2]))
        atrial2LJ.append(tools.midPt(atrialShort,
                                     atrialPts[i][2], atrialPts[i][3]))  #
        atrialJ4L.append(tools.midPt(atrialShort,
                                     atrialPts[i][3], atrialPts[i][4]))  #
        atrial2R4L.append(tools.midPt(atrialShort,
                                      atrialPts[i][4], atrialPts[i][5]))
        atrial3R2R.append(tools.midPt(atrialShort,
                                      atrialPts[i][5], atrialPts[i][0]))
        # Smooth
        atrial3Right.append(tools.midPt(atrialShort,
                                        atrial4R3R[-1], atrial3R2R[-1]))
        atrial4Right.append(tools.midPt(atrialShort,
                                        atrial2L4R[-1], atrial4R3R[-1]))
        atrial2Left.append(tools.midPt(atrialShort,
                                       atrial2LJ[-1], atrial2L4R[-1]))
        atrialJoint.append(tools.midPt(atrialShort,
                                       atrialJ4L[-1], atrial2LJ[-1]))
        atrial4Left.append(tools.midPt(atrialShort,
                                       atrialJ4L[-1], atrial2R4L[-1]))
        atrial2Right.append(tools.midPt(atrialShort,
                                        atrial2R4L[-1], atrial3R2R[-1]))
    if len(atrialJoint) <= 1:
        sys.exit('The short axis contours are too short')
    # Additional layer or cut
    # Two layers 8 mm thickness
    if np.dot(atrialJoint[-1]-atrialJoint[0], axisDir) < 16:
        dis = 16-np.dot(atrialJoint[-1]-atrialJoint[0], axisDir)
        atrial3Right.append(atrial3Right[-1]+dis*axisDir)
        atrial4Right.append(atrial4Right[-1]+dis*axisDir)
        atrial2Left.append(atrial2Left[-1]+dis*axisDir)
        atrialJoint.append(atrialJoint[-1]+dis*axisDir)
        atrial4Left.append(atrial4Left[-1]+dis*axisDir)
        atrial2Right.append(atrial2Right[-1]+dis*axisDir)
        atrial4R3R.append(atrial4R3R[-1]+dis*axisDir)
        atrial2L4R.append(atrial2L4R[-1]+dis*axisDir)
        atrial2LJ.append(atrial2LJ[-1]+dis*axisDir)
        atrialJ4L.append(atrialJ4L[-1]+dis*axisDir)
        atrial2R4L.append(atrial2R4L[-1]+dis*axisDir)
        atrial3R2R.append(atrial3R2R[-1]+dis*axisDir)

    ax.plot([row[0] for row in atrial3Right], [row[1] for row in atrial3Right],
            [row[2] for row in atrial3Right], '-x')
    ax.plot([row[0] for row in atrial4Right], [row[1] for row in atrial4Right],
            [row[2] for row in atrial4Right], '-x')
    ax.plot([row[0] for row in atrial2Left], [row[1] for row in atrial2Left],
            [row[2] for row in atrial2Left], '-x')
    ax.plot([row[0] for row in atrialJoint], [row[1] for row in atrialJoint],
            [row[2] for row in atrialJoint], '-x')
    ax.plot([row[0] for row in atrial4Left], [row[1] for row in atrial4Left],
            [row[2] for row in atrial4Left], '-x')
    ax.plot([row[0] for row in atrial2Right], [row[1] for row in atrial2Right],
            [row[2] for row in atrial2Right], '-x')
    ax.plot([row[0] for row in atrial4R3R], [row[1] for row in atrial4R3R],
            [row[2] for row in atrial4R3R], '-x')
    ax.plot([row[0] for row in atrial2L4R], [row[1] for row in atrial2L4R],
            [row[2] for row in atrial2L4R], '-x')
    ax.plot([row[0] for row in atrial2LJ], [row[1] for row in atrial2LJ],
            [row[2] for row in atrial2LJ], '-x')
    ax.plot([row[0] for row in atrialJ4L], [row[1] for row in atrialJ4L],
            [row[2] for row in atrialJ4L], '-x')
    ax.plot([row[0] for row in atrial2R4L], [row[1] for row in atrial2R4L],
            [row[2] for row in atrial2R4L], '-x')
    ax.plot([row[0] for row in atrial3R2R], [row[1] for row in atrial3R2R],
            [row[2] for row in atrial3R2R], '-x')

    refPt = atrialJoint[0]+8*axisDir
    numPTS = 10  # MUST BE EVEN

    iLA3Right = tools.trimLine(atrial3Right, refPt, axisDir, numPTS)
    iLA4R3R = tools.trimLine(atrial4R3R, refPt, axisDir, numPTS)
    iLA4Right = tools.trimLine(atrial4Right, refPt, axisDir, numPTS)
    iLA2L4R = tools.trimLine(atrial2L4R, refPt, axisDir, numPTS)
    iLA2Left = tools.trimLine(atrial2Left, refPt, axisDir, numPTS)
    iLA2LJ = tools.trimLine(atrial2LJ, refPt, axisDir, numPTS)
    iLAJoint = tools.trimLine(atrialJoint, refPt, axisDir, numPTS)
    iLAJ4L = tools.trimLine(atrialJ4L, refPt, axisDir, numPTS)
    iLA4Left = tools.trimLine(atrial4Left, refPt, axisDir, numPTS)
    iLA2R4L = tools.trimLine(atrial2R4L, refPt, axisDir, numPTS)
    iLA2Right = tools.trimLine(atrial2Right, refPt, axisDir, numPTS)
    iLA3R2R = tools.trimLine(atrial3R2R, refPt, axisDir, numPTS)
    ax.plot([row[0] for row in iLA3Right], [row[1] for row in iLA3Right],
            [row[2] for row in iLA3Right], '-')
    ax.plot([row[0] for row in iLA4Right], [row[1] for row in iLA4Right],
            [row[2] for row in iLA4Right], '-')
    ax.plot([row[0] for row in iLA2Left], [row[1] for row in iLA2Left],
            [row[2] for row in iLA2Left], '-x')
    ax.plot([row[0] for row in iLAJoint], [row[1] for row in iLAJoint],
            [row[2] for row in iLAJoint], '-x')
    ax.plot([row[0] for row in iLA4Left], [row[1] for row in iLA4Left],
            [row[2] for row in iLA4Left], '-x')
    ax.plot([row[0] for row in iLA2Right], [row[1] for row in iLA2Right],
            [row[2] for row in iLA2Right], '-x')
    ax.plot([row[0] for row in iLA4R3R], [row[1] for row in iLA4R3R],
            [row[2] for row in iLA4R3R], '-x')
    ax.plot([row[0] for row in iLA2L4R], [row[1] for row in iLA2L4R],
            [row[2] for row in iLA2L4R], '-x')
    ax.plot([row[0] for row in iLA2LJ], [row[1] for row in iLA2LJ],
            [row[2] for row in iLA2LJ], '-x')
    ax.plot([row[0] for row in iLAJ4L], [row[1] for row in iLAJ4L],
            [row[2] for row in iLAJ4L], '-x')
    ax.plot([row[0] for row in iLA2R4L], [row[1] for row in iLA2R4L],
            [row[2] for row in iLA2R4L], '-x')
    ax.plot([row[0] for row in iLA3R2R], [row[1] for row in iLA3R2R],
            [row[2] for row in iLA3R2R], '-x')

    # Aorta
    ooDir = np.array(joint)-np.array(iCh3Left[0])
    ooZDir = np.cross(ch3Normal, ooDir)/np.linalg.norm(np.cross(
        ch3Normal, ooDir))
    if np.dot(ooZDir, axisDir) > 0:
        aorticZAxis = ooZDir
    else:
        aorticZAxis = -1*ooZDir

    aoPts = [joint, atrial2LJ[0], iCh3L2L[0],
             iCh3Left[0], iCh4L3L[0], atrialJ4L[0]]
    aoCenter = np.mean(aoPts, axis=0)

    aoNormalAngle = []
    for i in range(0, len(aoPts)):
        oo = np.array(aoPts[i])-np.array(aoCenter)
        ooNormal = np.cross(aorticZAxis, oo)/np.linalg.norm(np.cross(
            aorticZAxis, oo))
        aoNormalAngle.append(ooNormal)

    ax.plot([row[0] for row in aoPts], [row[1] for row in aoPts],
            [row[2] for row in aoPts], '-<')

    ooMesh = [aoPts]
    for i in reversed(range(0, len(shorts[1]))):
        aorticShort = shorts[1][i]
        center = np.mean(aorticShort, axis=0)
        if aorticShort != [] and tools.l2s([apex, joint], center, axisDir) == []:
            '''
            ax.plot([row[0] for row in aorticShort],
                    [row[1] for row in aorticShort],
                    [row[2] for row in aorticShort], '-g')
            '''
            ax.scatter(center[0], center[1], center[2],
                       c='r', marker='o', s=50)
            interPts = []
            for n in range(0, len(aoPts)):
                oo = tools.l2sPt(aorticShort, center,
                                 aoNormalAngle[n], aorticZAxis)
                interPts.append(oo)
            ooMesh.append(interPts)

            ax.plot([row[0] for row in interPts], [row[1] for row in interPts],
                    [row[2] for row in interPts], '-r')

    # Add top layer

    oo = np.array(ooMesh[len(ooMesh)-1][0])-np.array(ooMesh[0][0])
    oo = np.dot(oo, aorticZAxis)

    if oo < 20:
        topPts = []
        refPt = (20-oo)*aorticZAxis+ooMesh[len(ooMesh)-1][0]
        for i in range(0, len(interPts)):
            ooPt = np.array((20-oo)*aorticZAxis)+np.array(interPts[i])
            topPts.append(ooPt)
        '''
        ax.plot([row[0] for row in topPts], [row[1] for row in topPts],
                [row[2] for row in topPts], '-xg')
        '''
        ooMesh.append(topPts)
    else:
        print 'Add'

    ao3Left = []
    ao4L3L = []
    aoJ4L = []
    aoJoint = []
    ao2LJ = []
    ao3L2L = []

    for i in range(0, len(ooMesh)):
        aoJoint.append(ooMesh[i][0])
        ao2LJ.append(ooMesh[i][1])
        ao3L2L.append(ooMesh[i][2])
        ao3Left.append(ooMesh[i][3])
        ao4L3L.append(ooMesh[i][4])
        aoJ4L.append(ooMesh[i][5])

    iAo3Left = tools.bSLine(ao3Left, len(ao3Left)-1, 50)
    iAo4L3L = tools.bSLine(ao4L3L, len(ao4L3L)-1, 50)
    iAoJ4L = tools.bSLine(aoJ4L, len(aoJ4L)-1, 50)
    iAoJoint = tools.bSLine(aoJoint, len(aoJoint)-1, 50)
    iAo2LJ = tools.bSLine(ao2LJ, len(ao2LJ)-1, 50)
    iAo3L2L = tools.bSLine(ao3L2L, len(ao3L2L)-1, 50)

    refPt = np.array(8*aorticZAxis)+np.array(joint)  # Cut
    ao3Left = tools.trimLine(iAo3Left, refPt, aorticZAxis, numPTS)
    ao4L3L = tools.trimLine(iAo4L3L, refPt, aorticZAxis, numPTS)
    aoJ4L = tools.trimLine(iAoJ4L, refPt, aorticZAxis, numPTS)
    aoJoint = tools.trimLine(iAoJoint, refPt, aorticZAxis, numPTS)
    ao2LJ = tools.trimLine(iAo2LJ, refPt, aorticZAxis, numPTS)
    ao3L2L = tools.trimLine(iAo3L2L, refPt, aorticZAxis, numPTS)

    ax.plot([row[0] for row in ao3Left], [row[1] for row in ao3Left],
            [row[2] for row in ao3Left], '-x')
    ax.plot([row[0] for row in ao4L3L], [row[1] for row in ao4L3L],
            [row[2] for row in ao4L3L], '-o')
    ax.plot([row[0] for row in aoJ4L], [row[1] for row in aoJ4L],
            [row[2] for row in aoJ4L], '-')
    ax.plot([row[0] for row in aoJoint], [row[1] for row in aoJoint],
            [row[2] for row in aoJoint], '-')
    ax.plot([row[0] for row in ao2LJ], [row[1] for row in ao2LJ],
            [row[2] for row in ao2LJ], '-')
    ax.plot([row[0] for row in ao3L2L], [row[1] for row in ao3L2L],
            [row[2] for row in ao3L2L], '-')
    # Avoid touching
    oo = tools.splitLines(aoJoint, iLAJoint)
    aoJoint = oo[0]
    iLAJoint = oo[1]
    oo = tools.splitLines(aoJ4L, iLAJ4L)
    aoJ4L = oo[0]
    iLAJ4L = oo[1]
    oo = tools.splitLines(ao2LJ, iLA2LJ)
    ao2LJ = oo[0]
    iLA2LJ = oo[1]
    
    venBone = [iCh2Left, iCh2Right, iCh3Left, iCh3Right, iCh4Left, iCh4Right,
               iCh3L2L, iCh2L4R, iCh4R3R, iCh3R2R, iCh2R4L, iCh4L3L]
    laBone = [iLA2Left, iLA2Right, iLAJoint, iLA3Right, iLA4Left, iLA4Right,
              iLA2LJ, iLA2L4R, iLA4R3R, iLA3R2R, iLA2R4L, iLAJ4L]
    aoBone = [ao3Left, ao4L3L, aoJ4L, aoJoint, ao2LJ, ao3L2L]
    return [venBone, laBone, aoBone, aorticZAxis]


def extractShortContours(fileName, totalNum, frameRate,
                         bSliceID, info, id, ax):
    venShorts = []
    aoShorts = []
    for i in range(0, totalNum/frameRate):
        shortName = "%s%04d.png" % (fileName, frameRate*i+id)
        colorPixel = p2dr1.loadContour(shortName, ['RED', 'LIME'], info[i+5])
        venPts = p2dr1.digitize3D(colorPixel[0], info[i+5])
        venShorts.append(venPts)
        aoPts = p2dr1.digitize3D(colorPixel[1], info[i+5])
        aoShorts.append(aoPts)

    venContours = []
    aoContours = []
    for i in range(0, bSliceID+1):
        venContour = tools.shortContour(venShorts[i], info[i+5], 100)
        venContours.append(venContour)
        '''
        if venContour != []:
            ax.plot([row[0] for row in venContour],
                    [row[1] for row in venContour],
                    [row[2] for row in venContour], '-g')
        '''
        aoContour = tools.shortContour(aoShorts[i], info[i+5], 100)
        aoContours.append(aoContour)
        '''
        if aoContour != []:
            ax.plot([row[0] for row in aoContour],
                    [row[1] for row in aoContour],
                    [row[2] for row in aoContour], '-r')
        '''
    return [venContours, aoContours]

# Create Folder if necessary
folderName = 'Simulation'
if not os.path.exists(folderName):
    os.makedirs(folderName)

if os.name == 'posix':
    folderName = 'Simulation/STL'
else:
    folderName = 'Simulation\\STL'
if not os.path.exists(folderName):
    os.makedirs(folderName)


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
'''
print 'LONG AXIS'
diaID=int(raw_input('The End of Diastole ID: '))
sysID=int(raw_input(' The End of Systole ID: '))
lfRate=int(raw_input(' The Total Number of long frames: ')
bSliceID=int(raw_input('   The bottom Slice ID: '))
#tSliceID=int(raw_input('      The top Slice ID: '))
ratio= float(raw_input('      Key in the ratio: '))

print 'SHORT AXIS'
fileName=     raw_input('    Key the file template name: ')
totalNum= int(raw_input('Key the total number of images: '))
frameRate=int(raw_input('            Key the frame rate: '))
'''
print 'Check Information on Line 1177'
bSliceID = 10
ratio = 0.65
lfRate = 40
fileName = 'short/short'
totalNum = 420
frameRate = 30

esID = int(raw_input('The end of systole ID: '))
esID = 12
atrialVolume = float(raw_input('The volume of atrium: '))
# Load long axis information
lCMatrix = []
for i in range(0, lfRate):
    if os.name == 'posix':
        centerName = 'longInfo/cLine_{}.txt'.format(i)
    else:
        centerName = 'longInfo\\cLine_{}.txt'.format(i)
    if not os.path.exists(centerName):
        sys.exit('cLine does not exist')
    else:
        longCenter = np.loadtxt(centerName)
        lCMatrix.append(longCenter)

# Interpolate to short
sCMatrix = []
for i in range(0, frameRate):
    sCMatrix.append([])
for i in range(0, len(longCenter)):
    pts = []
    for j in range(0, lfRate):
        pts.append(lCMatrix[j][i])
    ooPts = pts+pts+pts
    ooPts.append(pts[0])
    fullPts = tools.itpl(ooPts, frameRate*lfRate*3, 4, lfRate, 2*lfRate)
    for k in range(0, frameRate):
        # sCMatrix[k][i] = fullPts[k*lfRate]
        sCMatrix[k].append(fullPts[k*lfRate])

# Mitral valve
if os.name == 'posix':
    posteriorName = 'longInfo/fullPosterior.txt'
    anteriorName = 'longInfo/fullAnterior.txt'
else:
    posteriorName = 'longInfo\\fullPosterior.txt'
    anteriorName = 'longInfo\\fullAnterior.txt'

if os.path.exists(posteriorName) and os.path.exists(anteriorName):
    ch3pAngleLen = np.loadtxt(posteriorName)
    ch3aAngleLen = np.loadtxt(anteriorName)
    print 'Include Mitral'
else:
    ch3pAngleLen = []
    ch3aAngleLen = []
    print 'Exclude Mitral'

# Start building
sID = int(raw_input('The starting frameID: '))
eID = int(raw_input('  The ending frameID: '))
for id in range(sID, eID+1):
    print 'Frame ID: {}'.format(id)
    figure = plt.figure()
    ax = figure.gca(projection='3d')
    longCenter = sCMatrix[id][0:len(sCMatrix[id])-7]
    keyPts = sCMatrix[id][len(sCMatrix[id])-8:len(sCMatrix[id])]
    cDir = info[bSliceID+5][3:6]
    rDir = info[bSliceID+5][6:9]

    verticalDir = np.array(keyPts[1])-np.array(keyPts[0])
    if np.dot(verticalDir, np.cross(rDir, cDir)) > 0:
        axisDir = np.cross(rDir, cDir)
    else:
        axisDir = -1*np.cross(rDir, cDir)

    ax.plot([row[0] for row in longCenter],
            [row[1] for row in longCenter],
            [row[2] for row in longCenter],
            '-')
    ax.scatter(keyPts[0][0], keyPts[0][1],
               keyPts[0][2], c='r', marker='o', s=50)
    ax.scatter(keyPts[1][0], keyPts[1][1],
               keyPts[1][2], c='g', marker='o', s=50)
    ax.scatter(keyPts[2][0], keyPts[2][1],
               keyPts[2][2], c='r', marker='o', s=50)
    ax.scatter(keyPts[3][0], keyPts[3][1],
               keyPts[3][2], c='r', marker='o', s=50)
    ax.scatter(keyPts[4][0], keyPts[4][1],
               keyPts[4][2], c='r', marker='o', s=50)
    ax.scatter(keyPts[5][0], keyPts[5][1],
               keyPts[5][2], c='r', marker='o', s=50)
    ax.scatter(keyPts[6][0], keyPts[6][1],
               keyPts[6][2], c='r', marker='o', s=50)
    ax.scatter(keyPts[7][0], keyPts[7][1],
               keyPts[7][2], c='r', marker='o', s=50)
    shorts = extractShortContours(fileName, totalNum, frameRate,
                                  bSliceID, info, id, ax)

    bone = restoreGeo(axisDir, keyPts, longCenter, shorts, bSliceID, info, ax)
    if ch3pAngleLen != [] and ch3aAngleLen != []:
        mitralAngleLen = [ch3pAngleLen[id], ch3aAngleLen[id]]
        mitral = reCoverMitral(mitralAngleLen, keyPts, axisDir, info, ax)
    else:
        mitral = [[], []]
    atrialData = GEO(ratio, keyPts, bone, mitral, info, axisDir, esID, id, ax)
    outName = "short%04d.pdf" % id
    plt.savefig(outName, format='pdf')
    # plt.show()
    if atrialVolume != 0:
        reGenerateAtrium(atrialData, atrialVolume, axisDir, esID, id)       
