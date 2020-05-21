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
if __name__=="__main__":
    hand=bodySegment("Hand",{'proximalJointCentre':np.array([1,3,3]),'distalJointCentre':np.array([2,3,3])})
    DeLeva_BSIP_Table_Data=pd.read_table('DeLevaTable.txt',header=None,sep='\s+')
    DeLeva_BSIP_Table = pd.DataFrame(DeLeva_BSIP_Table_Data.values,
               index = ['Head','Trunk','UPT','MPT','LPT','UpperArm',
                        'ForeArm','Hand','Thigh','Shank','Foot'],
               columns=pd.MultiIndex.from_product([['Mass','CMPosition','Sagittal_r','Transverse_r','Longitudinal_r'],['Female','Male']]))
    hand.BSIP_Calculator(20,'Male',DeLeva_BSIP_Table)
    print(hand.I_Sagittal)