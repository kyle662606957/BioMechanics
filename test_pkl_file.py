import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import math
from collections import deque
import  re
import matplotlib as mpl


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
        self.timeLableKinematicsList.append(0)
        self.velocityMassCenterList=deque(maxlen=5)       
        self.accelerationMassCenterList =deque(maxlen=5)
        self.Tranformation_GCS2LCS_List=deque(maxlen=5)
        self.angularVelocityVectorList=deque(maxlen=5)
        self.angularAccelatrationList=deque(maxlen=5)
        self.jointLoadFlag=False
        
    def segmentLengthCalculator(self,proximalJointCentre=None,distalJointCentre=None):
        if proximalJointCentre==None:
            proximalJointCentre=self.proximalJointCentre
            distalJointCentre=self.distalJointCentre
        self.segmentlength=np.sqrt(np.sum((proximalJointCentre-distalJointCentre)**2))  
    def BSIP_Calculator(self,totalMass,Gender):
        nameIndex=re.sub("left|right", '', self.name)
        self.coefMassCenterPosition=self.DeLeva_BSIP_Table["CMPosition"][Gender][nameIndex]
        self.segmentMass=self.DeLeva_BSIP_Table["Mass"][Gender][nameIndex]*totalMass/100.0
        self.CMPosition=self.proximalJointCentre + \
                        self.coefMassCenterPosition*(self.distalJointCentre-
                        self.proximalJointCentre)/100.0
        self.I_Sagittal=self.segmentMass*(self.segmentlength*self.DeLeva_BSIP_Table["Sagittal_r"][Gender][nameIndex]/100.0)**2
        self.I_Transverse=self.segmentMass*(self.segmentlength*self.DeLeva_BSIP_Table["Transverse_r"][Gender][nameIndex]/100.0)**2
        self.I_Longitudinal=self.segmentMass*(self.segmentlength*self.DeLeva_BSIP_Table["Longitudinal_r"][Gender][nameIndex]/100.0)**2
        self.inertialMatrix=np.array([[self.I_Sagittal,0,0],[0,self.I_Transverse,0],[0,0,self.I_Longitudinal]])
    def insertChildren(self,listChildren):
        self.listChildren+=listChildren
    def proximalJointLoadCalculator(self):
        if self.jointLoadFlag==True:
            loadFromAllChildren=np.zeros([6])
            for childSegment in self.listChildren:
                loadFromAllChildren+=childSegment.proximalJointLoadCalculator()
            forcesJointFromAllChildren=loadFromAllChildren[0:3]
            momentJointFromAllChildren=loadFromAllChildren[3:]
            self.forcesProximal=(self.accelerationMassCenter-g_gravity)*self.segmentMass-forcesJointFromAllChildren
            self.momentsProximal=-np.cross(self.distalJointCentre-self.proximalJointCentre,forcesJointFromAllChildren) \
                                -np.cross(self.CMPosition-self.proximalJointCentre,self.segmentMass*(g_gravity-self.accelerationMassCenter))\
                                +np.transpose(self.Tranformation_GCS2LCS_List[-3])@self.inertialMatrix@self.Tranformation_GCS2LCS_List[-3]@self.angularAccelatration-momentJointFromAllChildren
            return np.concatenate((self.forcesProximal,self.momentsProximal))
    def UpdateKinematicInformation(self,markerData,timeLable):        
        self.timeLableKinematicsList.append(timeLable)
        self.proximalJointCentre=markerData["proximalJointCentre"]
        self.distalJointCentre=markerData["distalJointCentre"]
        self.CMPosition=self.proximalJointCentre + \
                        self.coefMassCenterPosition*(self.distalJointCentre-
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
        ## to calculate the tranformation matrix from the global coordinate system to local system
        ## the  are then calculated
        if len(markerData)==2:
            z_Axis_LCS=(self.distalJointCentre-self.proximalJointCentre)\
                /np.linalg.norm(self.distalJointCentre-self.proximalJointCentre)
            y_Axis_LCS=np.array([0,z_Axis_LCS[-1],z_Axis_LCS[-2]])\
                /np.linalg.norm(np.array([0,z_Axis_LCS[-1],z_Axis_LCS[-2]]))
            x_Axis_LCS=np.cross(y_Axis_LCS,z_Axis_LCS)
            self.Tranformation_GCS2LCS=np.vstack((x_Axis_LCS,y_Axis_LCS,z_Axis_LCS))
            self.Tranformation_GCS2LCS_List.append(self.Tranformation_GCS2LCS)
        ## calculate the joint angular velocity
        if len(self.Tranformation_GCS2LCS_List)>=3:
            R_Delta=np.matmul(self.Tranformation_GCS2LCS_List[-1],self.Tranformation_GCS2LCS_List[-3])
            angleValue=(R_Delta[0,0]+R_Delta[1,1]+R_Delta[2,2]-1)/2
            if abs(angleValue)>1:
                angleValue=angleValue/abs(angleValue)
            delta=math.acos(angleValue)
            omiga=delta/(self.timeLableKinematicsList[-1]-self.timeLableKinematicsList[-3])
            unitVector_V=np.array([R_Delta[1,2]-R_Delta[2,1],R_Delta[2,0]-R_Delta[0,2],R_Delta[0,1]-R_Delta[1,0]])\
                /(2.0*math.sin(delta))
            vector_V_GCS=np.matmul(np.transpose(self.Tranformation_GCS2LCS_List[-3]),unitVector_V)
            self.angularVelocityVector=omiga*vector_V_GCS
            self.angularVelocityVectorList.append(self.angularVelocityVector)
        if len(self.angularVelocityVectorList)>=3:
            self.jointLoadFlag=True
            self.angularAccelatration=(self.angularVelocityVectorList[-1]-self.angularVelocityVectorList[-3])\
                /(self.timeLableKinematicsList[-1]-self.timeLableKinematicsList[-3])
            self.angularAccelatrationList.append(self.angularAccelatration)
    def drawSegment(self,axToDraw):
        x=np.array([self.proximalJointCentre[0],self.distalJointCentre[0]])
        y=np.array([self.proximalJointCentre[1],self.distalJointCentre[1]])
        z=np.array([self.proximalJointCentre[2],self.distalJointCentre[2]])
        #axToDraw.redraw_in_frame(x,y,z)
        axToDraw.plot(x,y,z)

class humanBody:
    def __init__(self,totalBodyMass,gender,Visualization=True):
        self.totalBodyMass=totalBodyMass
        self.gender=gender
        if Visualization:
            self.fig=plt.figure()      
            self.axToDraw=self.fig.gca(projection='3d')
            
            
    def bodySegmentsGenerator(self,jointsCoordinatesDict,segmentDefinitionDict):
        self.segmentListDic={}
        for key in segmentDefinitionDict:
            exec('self.segmentListDic[\''+key+'\']=bodySegment(\''+key+'\',{\'distalJointCentre\':jointsCoordinatesDict[\''+segmentDefinitionDict[key][0]+'\'],\'proximalJointCentre\':jointsCoordinatesDict[\''+segmentDefinitionDict[key][1]+'\']})')   
        for segmentKey,segment in self.segmentListDic.items():
            segment.segmentLengthCalculator()
            segment.BSIP_Calculator(self.totalBodyMass,self.gender)
        self.segmentListDic['Trunk'].insertChildren([self.segmentListDic['Head'],self.segmentListDic['leftUpperArm'],self.segmentListDic['rightUpperArm']])
        self.segmentListDic['leftUpperArm'].insertChildren([self.segmentListDic['leftForeArm']])
        self.segmentListDic['rightUpperArm'].insertChildren([self.segmentListDic['rightForeArm']])  
    def bodySegmentsKinematicsUpdate(self,jointsCoordinatesDict,segmentDefinitionDict,timeLable):
        for segmentKey,segment in self.segmentListDic.items():
            segment.UpdateKinematicInformation({'proximalJointCentre':jointsCoordinatesDict[segmentDefinitionDict[segmentKey][1]],'distalJointCentre':jointsCoordinatesDict[segmentDefinitionDict[segmentKey][0]]},timeLable)
    def bodyLumbarLoadIntersegmental(self):
        return self.segmentListDic["Trunk"].proximalJointLoadCalculator()
    def drawBodySegment(self):
        self.axToDraw.cla()
        self.axToDraw.set_xlim([-1,1])
        self.axToDraw.set_ylim([-1,1])
        self.axToDraw.set_zlim([-1,1])
        self.axToDraw.set_xlabel("x")
        self.axToDraw.set_ylabel("y")
        self.axToDraw.set_zlabel("z")
        #self.fig.canvas.flush_events()
        for segmentKey,segment in self.segmentListDic.items():
            segment.drawSegment(self.axToDraw)        
        plt.pause(0.01)
        

results="/home/jindong/machineLearning/Human-Pose-Estimation/logs/eval_ibhgc_vol_softmax_VolumetricTriangulationNet@19.04.2021-17:33:12/checkpoints/0000/results.pkl"
with open(results, 'rb') as f:
    data = pickle.load(f)
#data keypoints_3d [frames,joints,positionxyz]
keypoints_position=data['keypoints_3d']

segmentDefinitionDict={'Head':['Head','ShoulderMid'],'Trunk':['ShoulderMid','HipCenter'],
                        'leftUpperArm':['L.Elbow','L.Shoulder'],'leftForeArm':['L.Wrist','L.Elbow'],#'leftHand':['L.Hand','L.Wrist'],'leftFoot':['L.Foot','L.Ankle'],
                        'rightUpperArm':['R.Elbow','R.Shoulder'],'rightForeArm':['R.Wrist','R.Elbow'],#'rightHand':['R.Hand','R.Wrist'],'rightFoot':['R.Foot','R.Ankle'],
                        'leftShank':['L.Ankle','L.Knee'],'leftThigh':['L.Knee','L.Hip'],
                        'rightShank':['R.Ankle','R.Knee'],'rightThigh':['R.Knee','R.Hip']}

listJoints=['R.Ankle','R.Knee','R.Hip','L.Hip','L.Knee','L.Ankle','HipCenter','Spine','ShoulderMid','Head',
            'R.Wrist','R.Elbow','R.Shoulder','L.Shoulder','L.Elbow','L.Wrist','nose']

jointsCoordinatesDict={}
initializationFlag=True
timeLable=0
body1=humanBody(75,'Male')
for frame in keypoints_position:
    timeLable+=1
    #jointCoordinates=np.array(list(map(float,line.split('\t'))))
    ##from the motion capture devices coodinates to the calculation coordinates
    #tranformationMatrix=np.array([[0,0,1],[-1,0,0],[0,1,0]])
    # for jointIndex,jointName in enumerate(listJoints):
    #     jointsCoordinatesDict[jointName]=tranformationMatrix@jointCoordinates[3*jointIndex:3*(jointIndex+1)]
    for jointIndex,jointName in enumerate(listJoints):
        jointsCoordinatesDict[jointName]=frame[jointIndex]         
    if initializationFlag:
        body1.bodySegmentsGenerator(jointsCoordinatesDict,segmentDefinitionDict)
        initializationFlag=False
    else:
        body1.bodySegmentsKinematicsUpdate(jointsCoordinatesDict,segmentDefinitionDict,timeLable)
        body1.drawBodySegment()
        load=body1.bodyLumbarLoadIntersegmental()
        try:
            print("Lumbar Load Calculation:\n  Forces (N):  {} \n  Moments(N.m):{} \n".format(load[0:3],load[3:]))
        except:
            pass
