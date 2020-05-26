import numpy as np
import pandas as pd
from collections import deque

g_gravity=np.array([0,0,-9.8])
class bodySegment:
    def __init__(self,name,markerData,Name_DeLeva_BSIP_Table='DeLevaTable.txt'):
        self.name=name
        self.proximalJointCentre=markerData["proximalJointCentre"]
        self.distalJointCentre=markerData["distalJointCentre"]
        self.segmentLengthCalculator()
        self.listChildren=[] 
        DeLeva_BSIP_Table_Data=pd.read_table(Name_DeLeva_BSIP_Table,header=None,sep='\s+')
        self.DeLeva_BSIP_Table = pd.DataFrame(DeLeva_BSIP_Table_Data.values,
               index = ['Head','Trunk','UPT','MPT','LPT','UpperArm',
                        'ForeArm','Hand','Thigh','Shank','Foot'],
               columns=pd.MultiIndex.from_product([['Mass','CMPosition','Sagittal_r','Transverse_r','Longitudinal_r'],['Female','Male']]))        
        self.CMPositionList=deque(maxlen=5)         
        self.timeLableKinematicsList=deque(maxlen=5)
        self.velocityMassCenterList=deque(maxlen=5)       
        self.accelerationMassCenterList =deque(maxlen=5)
    def segmentLengthCalculator(self,proximalJointCentre=None,distalJointCentre=None):
        if proximalJointCentre==None:
            proximalJointCentre=self.proximalJointCentre
            distalJointCentre=self.distalJointCentre
        self.segmentlength=np.sqrt(np.sum((proximalJointCentre-distalJointCentre)**2))  
    def BSIP_Calculator(self,totalMass,Gender):
        self.coefMassCenterPosition=self.DeLeva_BSIP_Table["CMPosition"][Gender][self.name]
        self.segmentMass=self.DeLeva_BSIP_Table["Mass"][Gender][self.name]*totalMass
        self.CMPosition=self.proximalJointCentre + \
                        self.coefMassCenterPosition*(self.proximalJointCentre-
                        self.proximalJointCentre)/100.0
        self.I_Sagittal=self.segmentMass*(self.segmentlength*self.DeLeva_BSIP_Table["Sagittal_r"][Gender][self.name]/100.0)**2
        self.I_Transverse=self.segmentMass*(self.segmentlength*self.DeLeva_BSIP_Table["Transverse_r"][Gender][self.name]/100.0)**2
        self.I_Longitudinal=self.segmentMass*(self.segmentlength*self.DeLeva_BSIP_Table["Longitudinal_r"][Gender][self.name]/100.0)**2
        self.inertialMatrix=np.array([[self.I_Sagittal,0,0],[0,self.I_Transverse,0],[0,0,self.I_Longitudinal]])
    def insertChildren(self,listChildren):
        self.listChildren+=listChildren
    def proximalJointLoadCalculator(self):
        loadFromAllChildren=0
        for childSegment in self.listChildren:
            loadFromAllChildren+=self.proximalJointLoadCalculator(childSegment)
        forcesJointFromAllChildren=loadFromAllChildren[0:3]
        momentJointFromAllChildren=loadFromAllChildren[3:]
        self.forcesProximal=(self.accelerationMassCenter-g_gravity)*self.segmentMass-forcesJointFromAllChildren
        self.momentsProximal=-np.cross(self.distalJointCentre-self.proximalJointCentre,forcesJointFromAllChildren) \
                              -np.cross(self.CMPosition-self.proximalJointCentre,self.segmentMass*(g_gravity-self.accelerationMassCenter))\
                             +np.matmul(self.inertialMatrix,self.angularVelocity)-momentJointFromAllChildren
        return np.concatenate(self.forcesProximal,self.momentsProximal)
    def UpdateKinematicInformation(self,markerData,timeLable):
        self.timeLableKinematicsList.append(timeLable)
        self.proximalJointCentre=markerData["proximalJointCentre"]
        self.distalJointCentre=markerData["distalJointCentre"]
        self.CMPosition=self.proximalJointCentre + \
                        self.coefMassCenterPosition*(self.proximalJointCentre-
                        self.proximalJointCentre)/100.0
        self.CMPositionList.append(self.CMPosition)        
        if len(self.CMPositionList)>=2:
            self.velocityCMPosition=(self.CMPositionList[-1]-self.CMPositionList[-2])\
                /(self.timeLableKinematicsList[-1]-self.timeLableKinematicsList[-2])
            self.velocityMassCenterList.append(self.velocityMassCenterList)
        if len(self.CMPositionList)>=4:
            self.accelerationMassCenter=(((self.CMPositionList[-1]-self.CMPositionList[-3])\
                /(self.timeLableKinematicsList[-1]-self.timeLableKinematicsList[-3]))-\
                    ((self.CMPositionList[-2]-self.CMPositionList[-4])\
                /(self.timeLableKinematicsList[-2]-self.timeLableKinematicsList[-4])))\
                    /(self.timeLableKinematicsList[-1]-self.timeLableKinematicsList[-2])
            self.accelerationMassCenterList.append(self.accelerationMassCenter)  
class humanBody:
    def __init__(self,totalBodyMass,gender):
        self.totalBodyMass=totalBodyMass
        self.gender=gender
    def bodySegmentsGenerator(self,jointsCoordinatesDict,segmentDefinitionDict):
        for key in segmentDefinitionDict:
            exec('self.'+key+'=bodySegment(\''+key+'\',{\'distalJointCentre\':jointsCoordinatesDict[\''+segmentDefinitionDict[key][0]+'\'],\'proximalJointCentre\':jointsCoordinatesDict[\''+segmentDefinitionDict[key][1]+'\']})')   


if __name__=="__main__":
    '''
    hand=bodySegment("Hand",{'proximalJointCentre':np.array([1,3,3]),'distalJointCentre':np.array([2,3,3])})
    DeLeva_BSIP_Table_Data=pd.read_table('DeLevaTable.txt',header=None,sep='\s+')
    DeLeva_BSIP_Table = pd.DataFrame(DeLeva_BSIP_Table_Data.values,
               index = ['Head','Trunk','UPT','MPT','LPT','UpperArm',
                        'ForeArm','Hand','Thigh','Shank','Foot'],
               columns=pd.MultiIndex.from_product([['Mass','CMPosition','Sagittal_r','Transverse_r','Longitudinal_r'],['Female','Male']]))
    hand.BSIP_Calculator(20,'Male',DeLeva_BSIP_Table)
    print(hand.I_Sagittal)
    '''
    with open('motionFile.trc','r') as f:
        dataLines=f.readlines()
        segmentDefinitionDict={'Head':['Head','ShoulderMid'],'Trunk':['ShoulderMid','HipCenter'],
                               'leftUpperArm':['L.Elbow','L.Shoulder'],'leftForeArm':['L.Wrist','L.Elbow'],'leftHand':['L.Hand','L.Wrist'],
                               'rightUpperArm':['R.Elbow','R.Shoulder'],'rightForeArm':['R.Wrist','R.Elbow'],'rightHand':['R.Hand','R.Wrist'],
                               'leftFoot':['L.Foot','L.Ankle'],'leftShank':['L.Ankle','L.Knee'],'leftThigh':['L.Knee','L.Hip'],
                               'rightFoot':['R.Foot','R.Ankle'],'rightShank':['R.Ankle','R.Knee'],'rightThigh':['R.Knee','R.Hip']}
        listJoints=['HipCenter','Spine','ShoulderMid','Head','L.Shoulder','L.Elbow','L.Wrist','L.Hand',
        'R.Shoulder','R.Elbow',	'R.Wrist','R.Hand','L.Hip','L.Knee','L.Ankle','L.Foot','R.Hip',
        'R.Knee','R.Ankle','R.Foot']
        jointsCoordinatesDict={}
        for line in dataLines:
            jointCoordinates=np.array(list(map(float,line.split('\t'))))
            for jointIndex,jointName in enumerate(listJoints):
                jointsCoordinatesDict[jointName]=jointCoordinates[3*jointIndex:3*(jointIndex+1)]
            body1=humanBody(75,'Male')
            body1.bodySegmentsGenerator(jointsCoordinatesDict,segmentDefinitionDict)
            print(body1)

                

