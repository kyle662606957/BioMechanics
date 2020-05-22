import numpy as np
import pandas as pd

class bodySegment:
    def __init__(self,name,markerData):
        self.name=name
        self.proximalJointCentre=markerData["proximalJointCentre"]
        self.distalJointCentre=markerData["distalJointCentre"]
        self.segmentLengthCalculator()
    def segmentLengthCalculator(self,proximalJointCentre=None,distalJointCentre=None):
        if proximalJointCentre==None:
            proximalJointCentre=self.proximalJointCentre
            distalJointCentre=self.distalJointCentre
        self.segmentlength=np.sqrt(np.sum((proximalJointCentre-distalJointCentre)**2))  
    def BSIP_Calculator(self,totalMass,Gender,DeLeva_BSIP_Table):
        self.segmentMass=DeLeva_BSIP_Table["Mass"][Gender][self.name]*totalMass
        self.CMPosition=self.proximalJointCentre + \
                        DeLeva_BSIP_Table["CMPosition"][Gender][self.name]*(self.proximalJointCentre-
                        self.proximalJointCentre)/100.0
        self.I_Sagittal=self.segmentMass*(self.segmentlength*DeLeva_BSIP_Table["Sagittal_r"][Gender][self.name]/100.0)**2
        self.I_Transverse=self.segmentMass*(self.segmentlength*DeLeva_BSIP_Table["Transverse_r"][Gender][self.name]/100.0)**2
        self.I_Longitudinal=self.segmentMass*(self.segmentlength*DeLeva_BSIP_Table["Longitudinal_r"][Gender][self.name]/100.0)**2
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

                

