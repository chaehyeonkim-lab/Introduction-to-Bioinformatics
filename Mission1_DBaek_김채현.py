#import sys, time
#sys.exit(sErrorMessage)
#print( time.ctime() )

import math

# file 처리
def file(file_name):
    
    handle = open(file_name, "r")
    handle.readline()
    fasta = handle.read().replace('\n','')
    handle.close()

    fasta = fasta.upper()
    
    return fasta

#nt counting class 만들기
class seq_count:

    def __init__(self, fasta, Nmer):
        self.fasta = fasta
        self.Nmer = int(Nmer)

    def counting(self, merDict):

        for i in range (len(self.fasta) - (self.Nmer-1)):
            nt = self.fasta[i:i + self.Nmer]
            
            if 'N' in nt:
                continue
            elif nt in merDict:
                merDict[nt] += 1
            else:
                merDict[nt] = 1
                
        return merDict
# end class

class ratio_cal:

    def __init__(self, merDict):
        self.merDict = merDict

    def ratio(self):
        total_count = 0
        ratio_sum = 0
        p_NT = {}

        for key in self.merDict:
            total_count += self.merDict[key]

        for key in self.merDict:
            freq = self.merDict[key] / total_count
            p_NT[key] = freq
            ratio_sum += freq

        return p_NT
#end class

def FC_Cal(twomer_p_NT, onemer_p_NT):
    cal_result = {}

    for key in twomer_p_NT:
        p_1stNT = key[0] #p_1stNT는 글자
        p_2ndNT = key[1]
        
        expected_freq = onemer_p_NT[p_1stNT] * onemer_p_NT[p_2ndNT]
        
        log2FC = math.log2(twomer_p_NT[key]/expected_freq)

        cal_result[key] = {"expected frequency": expected_freq, "log2FC": log2FC}

    return cal_result

def main(spieces_name, spieces, folder):
    onemer_Dict = {}
    twomer_Dict = {}
    
    for ch in spieces:
        file_name = '/Users/kimchaehyeon/Documents/seoul national univ./univ 4-2/생물학자를 위한 전산학 개론/'+folder+'/chr'+ch+'.fa'
        fasta = file(file_name)

        onemer = seq_count(fasta, 1)
        twomer = seq_count(fasta, 2)

        onemer_Dict = onemer.counting(onemer_Dict)
        twomer_Dict = twomer.counting(twomer_Dict)

    onemer_ratio = ratio_cal(onemer_Dict)
    twomer_ratio = ratio_cal(twomer_Dict)
    
    pVal_onemer = onemer_ratio.ratio() #1mer 비율 계산
    pVal_twomer = twomer_ratio.ratio() #2mer 비율 계산
        
    freq_cal = FC_Cal(pVal_twomer, pVal_onemer)

    print (spieces_name)
    for key in sorted(onemer_Dict):
        print(key, onemer_Dict[key], pVal_onemer[key], sep='\t')
    for key in sorted(twomer_Dict):
        print (key, twomer_Dict[key], pVal_twomer[key], freq_cal[key]["expected frequency"], freq_cal[key]["log2FC"], sep = '\t')            

C_elegans = ['I','II','III','IV','V','X']
D_melanogaster = ['2L','2R','3L','3R','X','Y']    
Chicken = [str(i) for i in range (1, 29)] + ['W','Z']
Human = [str(i) for i in range (1, 23)] + ['M','X','Y']

main("C_elegans", C_elegans, 'ce11ChrFiles')
#main(D_melanogaster, 'dm6ChrFiles')
#main(Chicken, 'galGal6ChrFiles')
#main(Human, 'hg38ChrFilesIncludingChrM')
